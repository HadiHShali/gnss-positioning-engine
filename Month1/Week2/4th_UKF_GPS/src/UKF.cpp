//==============================================================//
// Week 2, Day 4: Unecented Kalman Filter - GPS Altitude Fusion //
//==============================================================//

// No Jacobian - Sigma Points handle the Nonlinearity

// SECTION 0: Headers
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <Eigen/Dense>
#include <vector>
#include<iomanip>

using namespace std;
using namespace Eigen;

// SECTION 1: UKF Class
class UKF
{
public:
	int n;								// dimension of the state vector
	double alpha, beta, kappa, lambda;  // tuning parameters
	VectorXd x;							// state vector
	MatrixXd P;							// State Covariance Matrix
	MatrixXd Q;							// Process noise Matrix
	MatrixXd R;							// Measurement Covariance Matrix
	VectorXd Wm, Wc;					// Weights for mean and covariance

	// Constructor
	UKF(int state_dim)
		:n(state_dim)
	{
		// tuning parameters:
		alpha  = 1e-3;
		beta   = 2.0;
		kappa  = 0.0;
		lambda = alpha * alpha * (n + kappa) - n;

		// initialize the state 
		// Altitude = 1150.0 m, vertical velocity = 0.0
		x = VectorXd(n);  // create an empty vector of size n*1 (2*1)
		x << 1150.0, 0.0; // fill it: x(0)=1150.0, x(1)=0.0
		// Covariances 
		P = MatrixXd::Identity(n,n) * 500.0;  // High initial uncertainty in both states
		Q = MatrixXd::Identity(n, n); Q(0, 0) = 1.0; Q(1, 0) = 0.1; // Altitude drifts more than velocity per epoch
		R = MatrixXd::Identity(1, 1) * 225.0; // GPS altitude std = 15m -> 225m^2 variance

		// Weights for mean and variance
		Wm = VectorXd(2 * n + 1);
		Wm(0) = lambda / (n + lambda);
		Wc = VectorXd(2 * n + 1);
		Wc(0) = lambda / (n + lambda) + (1 - alpha * alpha + beta);
		for (int i = 1; i < 2 * n + 1; i++)
		{
			Wm(i) = 1.0 / (2 * (n + lambda));
			Wc(i) = 1.0 / (2 * (n + lambda));
		}
	}

	// Generate Sigma Points
	// σ₀=x̂, σᵢ=x̂+col_i(√((n+λ)P)), σᵢ₊ₙ=x̂-col_i(√((n+λ)P))
	MatrixXd sigmapoints()
	{
		MatrixXd sigma(n, 2 * n + 1); // an empty matrix of 2*5 
		sigma.col(0) = x; // the first column in going to be filled with the state vector
		// Cholesky Decomposition of scaled covariacne
		// LLT solver — Eigen's Cholesky method
		// matrixL() extract lower triangular L
		
		MatrixXd L = ((n + lambda) * P).llt().matrixL(); 
		for (int i = 0; i < n ; i++)
		{
			sigma.col(i + 1) = x + L.col(i); // positive spread
			sigma.col(i + n + 1) = x - L.col(i); // Negative spread
		}
		return sigma;
	}

	// Predict
	// Process model: altitude stays, velocity damps slightly
	// it describes how the state evolves from one time step to the next.
	
	// s : current sigma point (a 2×1 vector: [altitude, velocity])
	// dt : time step between epochs (typically 1.0 second for GPS) 
	VectorXd processModel(VectorXd s, double dt) 										 
	{
		VectorXd s_new(n);
		s_new(0) = s(0) + s(1) * dt; //altitude += vel * dt
		s_new(1) = s(1) * 0.99; // velocity damps slightly
		return s_new;
	}
	// prediction
	void predict(double dt = 1.0)
	{
		MatrixXd sigma = sigmapoints();
		MatrixXd sigma_pred(n, 2 * n + 1);

		// propagate each sigma point through process model
		for (int i = 0; i < 2 * n + 1; i++)
		{
			sigma_pred.col(i) = processModel(sigma.col(i), dt);
		}

		// Weighted Mean
		// x̂⁻ = Σ Wm_i * σᵢ⁻
		x = VectorXd::Zero(n);
		for (int i = 0; i < 2 * n + 1; i++)
		{
			x += Wm(i) * sigma_pred.col(i);
		}

		// Weighted covariance: 
		// P⁻ = Σ Wc_i*(σᵢ⁻-x̂⁻)(σᵢ⁻-x̂⁻)ᵀ + Q
		P = Q;
		for (int i = 0; i < 2 * n + 1; i++)
		{
			VectorXd d = sigma_pred.col(i) - x;
			P += Wc(i) * d * d.transpose();
		}
	}

	// UPDATE 
	// Measurement model: GPS observes altitude directly h(x)=x(0)
	VectorXd measurementModel(VectorXd s) 
	{
		VectorXd z(1); z(0) = s(0);  // GPS measures altitude
		return z;
	}


	void update(double gps_alt)
	{
		MatrixXd sigma = sigmapoints();
		VectorXd z_meas(1);
		z_meas(0) = gps_alt;

		// Propagate sigma points through measurement model
		MatrixXd Z(1, 2 * n + 1);
		for (int i = 0; i < 2 * n + 1; i++) Z.col(i) = measurementModel(sigma.col(i));

		// Predicted measurement: ẑ = Σ Wm_i * Zᵢ
		VectorXd z_pred = VectorXd::Zero(1);
		for (int i = 0; i < 2 * n + 1; i++) z_pred += Wm(i) * Z.col(i);

		// Innovation covariance: Pzz = Σ Wc_i*(Zᵢ-ẑ)(Zᵢ-ẑ)ᵀ + R
		MatrixXd Pzz = R;
		for (int i = 0; i < 2 * n + 1; i++) {
			VectorXd dz = Z.col(i) - z_pred;
			Pzz += Wc(i) * dz * dz.transpose();
		}

		// Cross-covariance: Pxz = Σ Wc_i*(σᵢ-x̂)(Zᵢ-ẑ)ᵀ
		MatrixXd Pxz = MatrixXd::Zero(n, 1);
		for (int i = 0; i < 2 * n + 1; i++) {
			VectorXd dx = sigma.col(i) - x;
			VectorXd dz = Z.col(i) - z_pred;
			Pxz += Wc(i) * dx * dz.transpose();
		}

		// Kalman gain: K = Pxz * Pzz⁻¹  (no Jacobian!)
		MatrixXd K = Pxz * Pzz.inverse();

		// Update state and covariance
		x = x + K * (z_meas - z_pred);
		P = P - K * Pzz * K.transpose();
	}
	double altitude() { return x(0); }
	double velocity() { return x(1); }

	
};

// SECTION 2: MAIN Program
int main() 
{
	UKF ukf(2);  // 2D state: [altitude, velocity]
	double gps[] = { 1205.3,1178.9,1198.2,1183.7,1204.1,
					1176.5,1195.8,1188.3,1201.9,1185.4 };
	cout << fixed << setprecision(2);
	cout << "=== UKF: GPS Altitude Fusion (Tehran 1191m) ===" << endl << endl;
	cout << "Epoch  GPS(m)    UKF Alt(m)  UKF Vel(m/s)" << endl;
	for (int i = 0; i < 10; i++) {
		ukf.predict();
		ukf.update(gps[i]);
		cout << "  " << i + 1 << "    " << gps[i]
			<< "    " << ukf.altitude()
			<< "    " << ukf.velocity() << endl;
	}
	cout << endl << "Final altitude: " << ukf.altitude() << " m  (true: 1191 m)" << endl;
	return 0;
}
