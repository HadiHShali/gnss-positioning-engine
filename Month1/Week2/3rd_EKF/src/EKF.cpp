// ====================================================================//
// Week 2, Day 3: Extended Kalman Filter (EKF) - GPS Range Measurement //
// ====================================================================//

// State: 2D receiver position = [x, y]
// Measurement: pseudorange -> rho = sqrt((Xs-x)^2 + (Ys-y)^2) 
// Jacobian (H) is recomputed every epoch at current estimate

// SECTION 0 : Headers
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include<iomanip>
using namespace std;
using namespace Eigen;

// SECTION 1: EKF class -----------------------------------------------------
class EKF_GPS
{
public:
	Vector2d x; // state vector [x_pos, y_pos] (2*1)
	Matrix2d P; // Covariance Matrix (2*2)
	Matrix2d F; // State transition -> identity: position is constant (2*2)
	Matrix2d Q; // Process Noise (2*2)
	double R; // Measurement Noise -> scalar: one range per epoch
	Matrix2d I; // Identity Matrix (2*2)

	// Constructor
	EKF_GPS(double x0, double y0)
	{
		x = Vector2d(x0, y0);
		P = Matrix2d::Identity() * 500.0; // high initial uncertainty
		F = Matrix2d::Identity(); // F = I: position stays constant
		Q = Matrix2d::Identity() * 1.0; // small process noise
		R = 100.0; // GPS range noise: 10 m of std = 100m^2 of Variance
		I = Matrix2d::Identity();

	}

	// Prediction Step ----------------------------------------------------
	// x̂⁻ = F·x̂     P⁻ = F·P·Fᵀ + Q
	void predict()
	{
		x = F * x;
		P = F * P * F.transpose() + Q;
	}

	// Non Linear Mesurement Function -------------------------------------
	// h(x) = sqrt((Xs-x)^2 + (Ys-y)^2)
	double h(double Xs, double Ys)
	{
		double dx = Xs - x(0);
		double dy = Ys - x(1);
		return sqrt(dx * dx + dy * dy);

	}

	//Jacobian ------------------------------------------------------------
	// H = [dh/dx, dh/dy] = [-(Xs-x)/rho, -(Ys-y)/rho]
	RowVector2d jacobian(double Xs, double Ys)
	{
		double dx = Xs - x(0);
		double dy = Ys - x(1);
		double rho = sqrt(dx * dx + dy * dy);
		return RowVector2d(-dx / rho, -dy / rho); // 1*2 jacobian
	}

	// Update -------------------------------------------------------------
	// y = z - h(x̂⁻)              (nonlinear innovation)
	// K = P⁻·Hᵀ·(H·P⁻·Hᵀ+R)⁻¹   (using Jacobian H)
	// x̂ = x̂⁻ + K·y
	// P = (I - K·H)·P⁻
	void update(double range_measure, double Xs, double Ys)
	{
		// compute Jacobian at current estimate
		RowVector2d H = jacobian(Xs, Ys);

		//Innovation: Measured range minus predicted range
		double y = range_measure - h(Xs, Ys);

		// Kalman Gain
		Vector2d K = P * H.transpose() / (H * P * H.transpose() + R);

		// update state and covariance
		x = x + K * y;
		P = (I - K * H) * P;
	}

	double posX() { return x(0); }
	double posY() { return x(1); }

};


// SECTION 2: Main Program
int main()
{
	// True position of receiver
	double true_x = 3352834.0, true_y = 4076956.0;

	// Initialize the EKF with rough initial guess (500m in x and 300m in y)
	EKF_GPS ekf(true_x - 500.0, true_y - 300.0);

	// Satellite ECEF position in 2D
	double Xs = 15600000.0, Ys = 7540000.0;

	// True range from recieverr to satellite
	double true_range = sqrt(pow(Xs - true_x, 2) + pow(Ys - true_y, 2));
	
	cout << "=== EKF: GPS Range Measurement ===" << endl;
	cout << "True position: (" << true_x << ", " << true_y << ")" << endl;
	cout << "True range to satellite: " << true_range / 1000.0 << " km" << endl << endl;
	cout << "Epoch  Range(m)    EKF_X        EKF_Y        Error(m)" << endl;
	cout << "-----  ----------  -----------  -----------  --------" << endl;

	cout << fixed << setprecision(2);
	for (int i = 0; i < 15; i++)
	{
		//predict
		ekf.predict();

		//simulate noisy range measurement (10m std)
		double noise = (rand() % 2000 - 1000) / 100.0; // -10 to 10 m
		double z = true_range + noise;

		//update
		ekf.update(z, Xs, Ys);

		// Make sure it looks exactly like this
		double err = sqrt(pow(ekf.posX() - true_x, 2) +
			pow(ekf.posY() - true_y, 2));

		cout << "  " << i + 1
			<< "    " << z
			<< "    " << ekf.posX()
			<< "    " << ekf.posY()
			<< "    " << err << endl;
	}

	return 0;


}