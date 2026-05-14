//===================================================================//
// Week 2, Day 2: 2D Kalmna Filter - GPS Position + Velocity Tracking//
//===================================================================//
// 
// State = [x, y, xdot, ydot] (4*1)
// Measurement = [xGPS, yGPS] (2*1)

// SECTION 0: Headers
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// SECTION 1: Kalman Filter Class
class KalmanFilter2D
{
public:
	MatrixXd x; // state = [x, y, xdot, ydot] (4*1)
	MatrixXd P; // covariance (4*4)
	MatrixXd F; // state trnsition (4*4)
	MatrixXd H; // Measurement model (2*4)
	MatrixXd Q; // Process noise (4*4)
	MatrixXd R; // Measurement noise (2*2)
	MatrixXd I; // Identity (4*4)

	// Constructor (Initialization)-----------------------------
	KalmanFilter2D(double dt)
	{
		// Initial state: at origin, moving 5m/s East, 3 m/s North
		x = Vector4d(0.0, 0.0, 5.0, 3.0);

		// Initial covariance: high uncertainty
		P = MatrixXd::Identity(4, 4) * 100.0; 

		// F: constant Velocity Model - drived from kinematic
		// x_{k+1} = x_k + dt * xdot_k
		// y_{k+1} = y_k + dt * ydot_k
		F = MatrixXd::Identity(4, 4);
		F(0, 2) = dt; 
		F(1, 3) = dt;

		// H: GPS observes x and y positions only (NOT velocity)
		H = MatrixXd::Zero(2, 4);
		H(0, 0) = 1.0;
		H(1, 1) = 1.0;

		// Q: Process noise  - velocity is more uncertain than position
		Q = MatrixXd::Zero(4, 4);
		Q(0, 0) = 0.1; // x position noise
		Q(1, 1) = 0.1; // y position noise
		Q(2, 2) = 1.0; // x velocity noise
		Q(3, 3) = 1.0; // y position noise

		// R: Measurement Noise - 8m std -> 64m^2 variance
		R = MatrixXd::Identity(2, 2) * 64.0; 

		// I: Identity MAtrix
		I = MatrixXd::Identity(4, 4);
	} // ----------------------------------------------------------------
	// Predic Step
	// Formula: x̂⁻ = F·x̂      P⁻ = F·P·Fᵀ + Q
	void predict()
	{
		x = F * x;
		P = F * P * F.transpose() + Q;
	} //-----------------------------------------------------------------

	// Update Step
	// Formula:
	// K = P⁻·Hᵀ·(H·P⁻·Hᵀ+R)⁻¹
    // x̂ = x̂⁻ + K·(z − H·x̂⁻)
    // P = (I − K·H)·P⁻
	void update(double xGPS, double yGPS)
	{
		// measurement vector z (2*1)
		Vector2d z(xGPS, yGPS); 

		// Kalman Gain (4*2)
		MatrixXd K = P * H.transpose() * (H * P * H.transpose() + R).inverse();

		// Innovation: Difference between GPS positions and the predicted positions
		VectorXd y = z - H * x;

		// Update the state
		x = x + K * y;

		// Update the covariance
		P = (I - K * H) * P;
	} //-----------------------------------------------------------------------

	// Getters: is a function that reads a private or internal value and returns it to you
	double posX()
	{ 
		return x(0);
	}  // returns element [0] = East position

	double posY()
	{ 
		return x(1);
	}  // returns element [1] = North position
	double velX()
	{ 
		return x(2);
	}  // returns element [2] = East velocity
	double velY()
	{ 
		return x(3);
	}  // returns element [3] = North velocity

};

// SECTION 2: Main Program
int main()
{
	// GPS at 1Hz -> dt = 1 sec : We are basically initializing the matrices
	KalmanFilter2D kf(1.0);

	// Simulated GPS readings: Vehicle is moving East at 5m/s, North at 3m/s
	// True Positions: (105, 203), (110, 206), (115, 209), ...
	double xGPS[] = { 113.2, 97.6, 118.4, 106.7, 122.1,
				  108.9, 125.8, 117.3, 131.9, 123.4 };
	double yGPS[] = { 204.8, 198.2, 211.3, 205.7, 216.4,
					  207.5, 219.8, 214.3, 224.9, 218.6 };

	//
	cout << "=== 2D Kalman Filters: GPS Position Tracker ====" << endl;
	cout << "Vehicle moving East at 5 m/s anmd North at 3m/s" << endl << endl;
	cout << "Epoch   GPS X   GPS Y   KF X   KF Y   Vel X   Vel Y" << endl;
	cout << "-----   -----   -----   ----   ----   -----   -----" << endl;

	for (int i = 0; i < 10; i++)
	{
		kf.predict();					 // step 1: prediction
		kf.update(xGPS[i], yGPS[i]);     // step 2: update

		cout << "   " << i + 1
			<< "   " << xGPS[i] << "   " << yGPS[i]
			<< "   " << kf.posX() << "   " << kf.posY()
			<< "   " << kf.velX() << "   " << kf.velY() << endl;

	}

	cout << endl << "Final Estimated velocity:" << endl;
	cout << " East Vel: xdot: " << kf.velX() << " m/s (true: 5 m/s)" << endl;
	cout << " North Vel: ydot: " << kf.velY() << " m/s (true: 3 m/s)" << endl;
	return 0;
}