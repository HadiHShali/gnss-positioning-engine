//
// Week 2, Day 1: 1D Kalman Filter - GPS Altitude smoother
//

// SECTION 0: Headers
#define _USE_MATH_DEFINES
#include <iostream>
#include <Eigen/Dense>
using namespace std;
using namespace Eigen;

//
// SECTION 1: Kalman Filter Class
class KalmanFilter1D
{
public:
	// State vector:
	MatrixXd x; // state estimate: [altitude] (1*1)

	// Covariance Matrix:
	MatrixXd P; // Covariance: uncertainty in altitude (1*1)

	// System matrices:
	MatrixXd F; // state trnsition: F=[1] altitude is constant
	MatrixXd H; // measurement: H=[1] GPS observes altitude directly
	MatrixXd Q; // Process noise: Q=[1] small random drift
	MatrixXd R; // Measurement noise: R=[225] GPS std deviation in altitude is 15m
	MatrixXd I; // Identity Matrix: I=[1] 

	// Constructor (initialization)
	KalmanFilter1D(double x0, double P0, double q, double r) {
		x = MatrixXd(1, 1); x(0, 0) = x0; // initial altitude estimate
		P = MatrixXd(1, 1); P(0, 0) = P0; // initial uncertaity
		F = MatrixXd(1, 1); F(0, 0) = 1.0; // F=[1]: altitude stays constant
		H = MatrixXd(1, 1); H(0, 0) = 1.0; // H=[1]: GPS observes altitude directly
		Q = MatrixXd(1, 1); Q(0, 0) = q; // Process noise (1m^2 per epoch)
		R = MatrixXd(1, 1); R(0, 0) = r; // measurement noise (15m std = 225 var)
		I= MatrixXd::Identity(1, 1);  // I = [1] for 1D
	}

	// Predic Step
	// Formula: x̂⁻ = F·x̂     P⁻ = F·P·Fᵀ + Q
	void predict()
	{
		x = F * x;  // for altitude: x stays the same
		P = F * P * F.transpose() + Q;  // uncertainty grows by Q=1 each epoch. 
	}

	// Update step
	// Equations:  K = P⁻·Hᵀ·(H·P⁻·Hᵀ+R)⁻¹
	//             x̂ = x̂⁻ + K·(z − H·x̂⁻)
	//             P = (I − K·H)·P⁻
	void update(double measurement)
	{
		MatrixXd z(1, 1); 
		z(0, 0) = measurement;

		// Kalman Gain: how much to trust GPS vs Model
		MatrixXd K = P * H.transpose() * (H * P * H.transpose() + R).inverse();

		// Innovation: difference between GPS readings vs predicted altitude
		MatrixXd y = z - H * x;

		// Correct state
		x = x + K * y;
		
		// Update uncertainty
		P = (I - K * H) * P; 

	}
	double altitude()
	{
		return x(0, 0); 
	}
	double uncertainty()
	{
		return sqrt(P(0, 0));
	}

};


// SECTION 2: Main Program
int main()
{
	// initialize the KF: x0 = 1150m 
	//					  P0 = 500
	//				      Q = 1
	//                    R = 225
	KalmanFilter1D kf(1150.0, 500.0, 1.0, 225.0);

	// Simulated GPS altitude readings for Tehran (true = 1191m)
	double gps[] = { 1205.3, 1178.9, 1198.2, 1183.7, 1204.1,
					1176.5, 1195.8, 1188.3, 1201.9, 1185.4 };

	cout << "=== 1D Kalman Filter: GPS Altitude (Tehran) ===" << endl;
	cout << "True altitude: 1191.0 m" << endl << endl;
	cout << "Epoch   GPS(m)   KF Est(m)   Uncertaity" << endl; 
	cout << "-----   ------   ---------   ----------" << endl;
	for (int i = 0; i < 10; i++)
	{
		kf.predict(); // step1: propagate forward
		kf.update(gps[i]); // step2: correct with measurement

		cout << "   " << i + 1 << "    " << gps[i] << "    " << kf.altitude() << "    " << kf.uncertainty() << " m" << endl;
	}

	cout << endl;
	cout << "Final KF estimate : " << kf.altitude() << " m" << endl;
	cout << "True Altitude: 1191.0m " << endl;
	cout << "Final Error: " << abs(kf.altitude() - 1191.0) << " m" << endl;
	cout << "GPS Noise (std): 15.0m -> KF reduces to ~" << kf.uncertainty() << "m" << endl;
	return 0;
}