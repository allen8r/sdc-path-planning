#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

vector<double> shift_reference(double x, double y, double ref_x, double ref_y, double ref_yaw) {
	// Shift car reference angle to 0 degrees
	double shift_x = x - ref_x;
	double shift_y = y - ref_y;

	shift_x = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
	shift_y = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));

	return {shift_x, shift_y};
}


/**
 * Returns a spline s. A vector of doubles for ref_x, ref_y, and ref_yaw is passed back to the &ref_states parameter
 * 
 */
tk::spline get_trajectory_spline(const double car_x, const double car_y, const double car_yaw,
																 const int prev_size, const vector<double> previous_path_x,
																 const vector<double> previous_path_y, const vector<double> &map_waypoints_s,
																 const vector<double> &map_waypoints_x, const vector<double> &map_waypoints_y,
																 const double S_BFFR, const double LN_WTH, const double LN_HLF,
																 const double car_s, const int &lane, vector<double> &ref_states) {
	
	cout << "get_trajectory_spline()..." << endl;
	vector<double> ptsx;
	vector<double> ptsy;

	// reference x, y, and yaw states
	// ref will be either the starting point of the current car position state or the preivous path's end point
	double ref_x = car_x;
	double ref_y = car_y;
	double ref_yaw = deg2rad(car_yaw);

	// if prev size is almost 0, then use the car's state as the starting reference
	if (prev_size < 2)
	{
		// Use 2 points that make the path tangent to the car
		double prev_car_x = car_x - cos(car_yaw);
		double prev_car_y = car_y - sin(car_yaw);

		ptsx.push_back(prev_car_x);
		ptsx.push_back(car_x);

		ptsy.push_back(prev_car_y);
		ptsy.push_back(car_y);
	}
	// Use the prev path's end point as the starting ref
	else
	{
		// Redefine reference state as the prev path's end point
		ref_x = previous_path_x[prev_size-1];
		ref_y = previous_path_y[prev_size-1];

		double ref_x_prev = previous_path_x[prev_size-2];
		double ref_y_prev = previous_path_y[prev_size-2];
		ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

		// Use 2 points that make the path tangent to the prev path's end point
		ptsx.push_back(ref_x_prev);
		ptsx.push_back(ref_x);

		ptsy.push_back(ref_y_prev);
		ptsy.push_back(ref_y);
	}

	ref_states.push_back(ref_x);
	ref_states.push_back(ref_y);
	ref_states.push_back(ref_yaw);

	// In Frenet add evenly 30m spaced points ahead of the starting reference
	double s_inc = S_BFFR;
	vector<double> next_wp0 = getXY(car_s+s_inc*1, (LN_HLF+LN_WTH*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
	vector<double> next_wp1 = getXY(car_s+s_inc*2, (LN_HLF+LN_WTH*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
	vector<double> next_wp2 = getXY(car_s+s_inc*3, (LN_HLF+LN_WTH*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

	ptsx.push_back(next_wp0[0]);
	ptsx.push_back(next_wp1[0]);
	ptsx.push_back(next_wp2[0]);

	ptsy.push_back(next_wp0[1]);
	ptsy.push_back(next_wp1[1]);
	ptsy.push_back(next_wp2[1]);

	// Transform points by shifting frame of reference to make math work easier
	for (int i = 0; i < ptsx.size(); i++)
	{
		// Shift car reference angle to 0 degrees
		double shift_x = ptsx[i] - ref_x;
		double shift_y = ptsy[i] - ref_y;

		ptsx[i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
		ptsy[i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
	}

	// Create a spline
	tk::spline s;
	
	// Set (x, y) points to the spline
	s.set_points(ptsx, ptsy);

	return s;
	
}

/**
 * Returns a trajectory of points. The trajectory is stored in a vector of vector of doubles. The first
 * vector of doubles are the x-values of the trajectory points, while the second vector of doubles is the
 * vector of y-values of the trajectory points.
 * 
 */
vector<vector<double>> gen_trajectory(const int &lane, const double car_x, const double car_y, 
																			const double car_yaw, const double car_s, vector<double> previous_path_x,
																			vector<double> previous_path_y, const int prev_size, double &ref_vel,
																			const vector<double> &map_waypoints_s, const vector<double> &map_waypoints_x,
  																		const vector<double> &map_waypoints_y, const double S_BFFR,
																			const double LN_WTH, const double LN_HLF, const double TIME_SLC) {
	// Smooth out the trajectory curve
	cout << "gen_trajectory(): lane: " << lane << endl;
	vector<double> ref_states;

	tk::spline s = get_trajectory_spline(car_x, car_y, car_yaw, prev_size, previous_path_x,
																			previous_path_y, map_waypoints_s, map_waypoints_x, map_waypoints_y,
																			S_BFFR, LN_WTH, LN_HLF, car_s, lane, ref_states);
	double ref_x = ref_states[0];
	double ref_y = ref_states[1];
	double ref_yaw = ref_states[2];

	// Define the actual (x, y) points for use in the path planner (points for the *future* path)
	vector<double> next_x_vals;
	vector<double> next_y_vals;

	// Start with all of the prev path points from the last cycle
	for (int i = 0; i < previous_path_x.size(); i++)
	{
		next_x_vals.push_back(previous_path_x[i]);
		next_y_vals.push_back(previous_path_y[i]);
	}

	// Calculate the spline points that will make the car travel at the desired velocity (Add triangle visual)
	double target_x = S_BFFR;
	double target_y = s(target_x);
	double target_dist = sqrt((target_x * target_x) + (target_y * target_y));

	double x_add_on = 0;

	// Fill up the rest of the path planner after filling it with prev points; always output 50 points here
	for (int i = 0; i <= 50-previous_path_x.size(); i++)
	{
		double N = (target_dist / (TIME_SLC * ref_vel/2.24)); // TODO: make magic numbers into constants (2.24 m/s)
		double x_point = x_add_on + (target_x / N);
		double y_point = s(x_point);

		x_add_on = x_point;

		double x_ref = x_point;
		double y_ref = y_point;

		// Rotate the point back to normal frame of reference after earlier shift
		x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
		y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

		x_point += ref_x;
		y_point += ref_y;

		next_x_vals.push_back(x_point);
		next_y_vals.push_back(y_point);
	}

	vector<vector<double>> trajectory;
	trajectory.push_back(next_x_vals);
	trajectory.push_back(next_y_vals);
	
	return trajectory;
}

bool safe_to_change_lane(const int target_lane, const vector<vector<double>> sensor_fusion,
													const double LN_WTH, const double prev_size,
													const double TIME_SLC, const double car_s, const double LC_BFFR) {
	bool safe_to_go = true;

	for (int j = 0; j < sensor_fusion.size(); j++)
	{
		float d_pos = sensor_fusion[j][6];
		
		if ((LN_WTH*target_lane) < d_pos && d_pos < (LN_WTH*(target_lane + 1)))
		{
			double vx = sensor_fusion[j][3];
			double vy = sensor_fusion[j][4];
			double check_speed = sqrt((vx * vx) + (vy * vy)); // find the magnitude of the velocity
			double check_car_s = sensor_fusion[j][5];

			// *Using prev points, project s value outward in time; look into future*
			check_car_s += ((double)prev_size * TIME_SLC * check_speed);
			
			/*
			cout << "---------------------------------------" << endl;
			cout << "                 d_pos: " << d_pos << endl;
			cout << "                 car_s: " << car_s << endl;
			cout << "           check_car_s: " << check_car_s << endl;
			cout << " (check_car_s - car_s): " << (check_car_s - car_s) << endl;
			cout << " (car_s - check_car_s): " << (car_s - check_car_s) << endl;
			cout << "---------------------------------------" << endl;
			*/

			if (((check_car_s > car_s) && ((check_car_s - car_s) < LC_BFFR)) ||
				((check_car_s < car_s) && ((car_s - check_car_s) < LC_BFFR)))
			{
				cout << "Car in the way..." << endl;
				safe_to_go = false;
			}
		}
	}
	return safe_to_go;
}

int change_lane(const int lane, const vector<vector<double>> sensor_fusion,
								const double LN_WTH, const int prev_size, const double TIME_SLC,
								const double car_s, const double LC_BFFR, const int LN_CTR,
								const int LN_LFT, const int LN_RGT) {
	
	int new_lane = lane;
  int target_lane;
	
	if (lane == LN_CTR)
	{
		// Check left lane
		if (safe_to_change_lane(lane-1, sensor_fusion, LN_WTH,
														prev_size, TIME_SLC, car_s, LC_BFFR))
		{
			new_lane = lane - 1;
		}
		// Check right lane
		else if (safe_to_change_lane(lane+1, sensor_fusion, LN_WTH,
														prev_size, TIME_SLC, car_s, LC_BFFR))
		{
			new_lane = lane + 1;
		}

	}
	else if (lane == LN_RGT)
	{
		// change one lane over to LEFT
		if (safe_to_change_lane(lane-1, sensor_fusion, LN_WTH,
														prev_size, TIME_SLC, car_s, LC_BFFR))
		{
			new_lane = lane - 1;
		}
	}
	else if (lane == LN_LFT)
	{
		// change one lane over to RIGHT
		if (safe_to_change_lane(lane+1, sensor_fusion, LN_WTH,
														prev_size, TIME_SLC, car_s, LC_BFFR))
		{
			new_lane = lane + 1;		
		}
	}

  cout << "Change lane to: " << new_lane << endl;
	return new_lane;
}


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

	double ref_vel = 0.0; // start at zero
	int lane = 1; // middle lane

  h.onMessage([&ref_vel, &lane, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;

		// Lane designations
		const int LN_LFT = 0;
		const int LN_CTR = 1;
		const int LN_RGT = 2;
		
		// Lane 'd' dimensions
		const double LN_WTH = 4;
		const double LN_HLF = LN_WTH / 2;
		
		const double VEL_MAX = 49.5;
		const double S_BFFR = 30.0; // 30 m
		//const double ACC_INC = 0.224; // acceleration increment, ~0.5 mph
		//const double ACC_INC = 0.447; // acceleration increment, ~1.0 mph
		const double ACC_INC = 0.5588; // acceleration increment, ~1.25 mph
		//const double ACC_INC = 0.67056; // acceleration increment, ~1.5 mph
		const double TIME_SLC = 0.02; // time slice in seconds (200 ms)
		const double LC_BFFR = 6.5; // 10 m


    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner (left-over points that were not used in the prev iteration in the simulator)
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

						/*
						for (int i = 0; i < sensor_fusion.size(); i++) {
								cout << "=====" << endl;
								cout << "id: " << sensor_fusion[i][0] << endl;
								cout << " x: " << sensor_fusion[i][1] << endl;
								cout << " y: " << sensor_fusion[i][2] << endl;
								cout << "vx: " << sensor_fusion[i][3] << endl;
								cout << "vy: " << sensor_fusion[i][4] << endl;
								cout << " s: " << sensor_fusion[i][5] << endl;
								cout << " d: " << sensor_fusion[i][6] << endl;
						}
						*/

          	json msgJson;
						int prev_size = previous_path_x.size();
						
						

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds


						// Handle sensor fusion (other cars' data)
						if (prev_size > 0)
						{
							car_s = end_path_s;
						}

						bool too_close = false;

						// Find ref_v
						for (int i= 0; i < sensor_fusion.size(); i++)
						{
							float d = sensor_fusion[i][6]; // d value tells which lane a car is in
							//if (d > (LN_HLF + LN_WTH*lane - LN_HLF) && d < (LN_HLF + LN_WTH*lane + LN_HLF))
							if (d > (LN_WTH*lane) && d < (LN_WTH*(lane + 1)))
							{
								double vx = sensor_fusion[i][3];
								double vy = sensor_fusion[i][4];
								double check_speed = sqrt((vx * vx) + (vy * vy)); // find the magnitude of the velocity
								double check_car_s = sensor_fusion[i][5];

								// *Using prev points, project s value outward in time; look into future*
								check_car_s += ((double)prev_size * TIME_SLC * check_speed); 

								// Check s values greater than ego car's and other car's s gap
								if ((check_car_s > car_s) && ((check_car_s - car_s) < S_BFFR)) // closer than 30m
								{
									// React because ego car is now too close to other car in the same lane.
									// Lower ref velocity so ego car avoids hitting car in front or
									// flag to change lane.
									//ref_vel = 29.5; // mph
									too_close = true;
									
									lane = change_lane(lane, sensor_fusion, LN_WTH, prev_size, TIME_SLC,
																			car_s, LC_BFFR, LN_CTR, LN_LFT, LN_RGT);
								}

							}
						}

						// Reduce/eliminate jerk by slowing and accelerating by increments up to
						// max velocity
						if (too_close)
						{
							ref_vel -= ACC_INC;
						}
						else if (ref_vel < VEL_MAX)
						{
							ref_vel += ACC_INC;
						}

						// BEGIN generating trajectory ============================================================

						vector<double> ptsx;
						vector<double> ptsy;

						// reference x, y, and yaw states
						// ref will be either the starting point of the current car position state or the preivous path's end point
						double ref_x = car_x;
						double ref_y = car_y;
						double ref_yaw = deg2rad(car_yaw);

						// if prev size is almost 0, then use the car's state as the starting reference
						if (prev_size < 2)
						{
							// Use 2 points that make the path tangent to the car
							double prev_car_x = car_x - cos(car_yaw);
							double prev_car_y = car_y - sin(car_yaw);

							ptsx.push_back(prev_car_x);
							ptsx.push_back(car_x);

							ptsy.push_back(prev_car_y);
							ptsy.push_back(car_y);
						}
						// Use the prev path's end point as the starting ref
						else
						{
							// Redefine reference state as the prev path's end point
							ref_x = previous_path_x[prev_size-1];
							ref_y = previous_path_y[prev_size-1];

							double ref_x_prev = previous_path_x[prev_size-2];
							double ref_y_prev = previous_path_y[prev_size-2];
							ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

							// Use 2 points that make the path tangent to the prev path's end point
							ptsx.push_back(ref_x_prev);
							ptsx.push_back(ref_x);

							ptsy.push_back(ref_y_prev);
							ptsy.push_back(ref_y);
						}

						// In Frenet add evenly 30m spaced points ahead of the starting reference
						double s_inc = S_BFFR;
						vector<double> next_wp0 = getXY(car_s+s_inc*1, (LN_HLF+LN_WTH*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
						vector<double> next_wp1 = getXY(car_s+s_inc*2, (LN_HLF+LN_WTH*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
						vector<double> next_wp2 = getXY(car_s+s_inc*3, (LN_HLF+LN_WTH*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

						ptsx.push_back(next_wp0[0]);
						ptsx.push_back(next_wp1[0]);
						ptsx.push_back(next_wp2[0]);

						ptsy.push_back(next_wp0[1]);
						ptsy.push_back(next_wp1[1]);
						ptsy.push_back(next_wp2[1]);

						// Transform points by shifting frame of reference to make math work easier
						for (int i = 0; i < ptsx.size(); i++)
						{
							// Shift car reference angle to 0 degrees
							double shift_x = ptsx[i] - ref_x;
							double shift_y = ptsy[i] - ref_y;

							ptsx[i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
							ptsy[i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
						}

						// Create a spline
						tk::spline s;
						
						// Set (x, y) points to the spline
						s.set_points(ptsx, ptsy);

						// Define the actual (x, y) points for use in the path planner (points for the *future* path)
						vector<double> next_x_vals;
						vector<double> next_y_vals;

						// Start with all of the prev path points from the last cycle
						for (int i = 0; i < previous_path_x.size(); i++)
						{
							next_x_vals.push_back(previous_path_x[i]);
							next_y_vals.push_back(previous_path_y[i]);
						}

						// Calculate the spline points that will make the car travel at the desired velocity (Add triangle visual)
						double target_x = S_BFFR;
						double target_y = s(target_x);
						double target_dist = sqrt((target_x * target_x) + (target_y * target_y));

						double x_add_on = 0;

						// Fill up the rest of the path planner after filling it with prev points; always output 50 points here
						for (int i = 0; i <= 50-previous_path_x.size(); i++)
						{
							double N = (target_dist / (TIME_SLC * ref_vel/2.24)); // TODO: make magic numbers into constants (2.24 m/s)
							double x_point = x_add_on + (target_x / N);
							double y_point = s(x_point);

							x_add_on = x_point;

							double x_ref = x_point;
							double y_ref = y_point;

							// Rotate the point back to normal frame of reference after earlier shift
							x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
							y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

							x_point += ref_x;
							y_point += ref_y;

							next_x_vals.push_back(x_point);
							next_y_vals.push_back(y_point);
						}
						// END generating trajectory ==============================================================

          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
