// kf_spp_processor.cpp: the executable that pipes trajectory.csv through your Kalman filter,
// RTS smoother, and Friedland bias estimator -> exports a 4-way comparison CSV

// ----------------------Section 0: Header + Setup------------------------------------------//
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "kalman_filter_gps.h"
#include "trajectory_reader.h"
#include "utils.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::Matrix3d;


int main() {
	cout << "=== KF SPP PROCESSOR ===" << endl;

	// Load raw SPP trajectory
	string in_csv = "../../data/trajectory_crctd.csv";
	auto epochs = readTrajectoryCsv(in_csv);
	cout << "Loaded " << epochs.size() << " raw epochs" << endl;

	if (epochs.size() < 4) {
		cerr << "Not enough epochs!" << endl;
		return 1;
	}


	// ------ initialize the KF with the first epoch-----------------------//
	//
	// static station -> small q_accel
	KalmanFilterGps kf(0.01);

	// Initial State from FIRST epoch (should be in ECEF)
	double x0, y0, z0;
	llaToEcef(epochs[0].lat_deg, epochs[0].lon_deg, epochs[0].height_m,
		x0, y0, z0);

	VectorXd x_init = VectorXd::Zero(8);
	x_init(0) = x0;   // position ECEF
	x_init(1) = y0;
	x_init(2) = z0;

	// Velocities and Clock terms all start at 0

	// Initialize the Weight Matrix
	MatrixXd P_init = MatrixXd::Identity(8, 8);
	P_init.block(0, 0, 3, 3) *= (50.0 * 50.0);   // 50m 1-sigma pos uncertainty
	P_init.block(3, 3, 3, 3) *= (1.0 * 1.0);      // 1 m/s velocity uncertainty
	P_init(6, 6) = 100.0;                          // clock bias uncertainty
	P_init(7, 7) = 1.0;                             // clock drift uncertainty


	kf.initialize(x_init, P_init);
	cout << "KF initialized at first epoch" << endl;


	double prev_t = epochs[0].t_gps;


	// Storage for RTS smoother
	vector<VectorXd> x_filt_arr, x_pred_arr;
	vector<MatrixXd> P_filt_arr, P_pred_arr, F_arr;

	x_filt_arr.reserve(epochs.size());
	x_pred_arr.reserve(epochs.size());
	P_filt_arr.reserve(epochs.size());
	P_pred_arr.reserve(epochs.size());
	F_arr.reserve(epochs.size());

	// Bias Estimator (Friedland)
	BiasEstimator biasEst;
	vector<Vector3d> bhat_history;
	vector<VectorXd> delta_history;   // <-- ADDED: per-epoch delta, needed for Variant B

	bhat_history.reserve(epochs.size());
	delta_history.reserve(epochs.size());

	// ------ Main loop-----------------------//
	//
	// Process Each Epoch Through Kalman Filter
	for (size_t i = 0; i < epochs.size(); i++)
	{
		const auto& e = epochs[i];

		// Step1: Convert lat/lon/h to ECEF (measurement z)
		double x, y, z;
		llaToEcef(e.lat_deg, e.lon_deg, e.height_m, x, y, z);

		Vector3d z_meas;
		z_meas(0) = x;
		z_meas(1) = y;
		z_meas(2) = z;


		// Step2: Predict (skip the first epoch)
		double dt = 15.0; // default for epoch 0
		if (i > 0) {
			dt = e.t_gps - prev_t;
			if (dt <= 0.0 || dt > 60.0) dt = 15.0;   // safety
			kf.predict(dt);
		}


		// SAVE prediction
		x_pred_arr.push_back(kf.getState());
		P_pred_arr.push_back(kf.getCovariance());


		// Reconstruct F used this step (matches buildF() in your class)
		MatrixXd F = MatrixXd::Identity(8, 8);
		F(0, 3) = dt; F(1, 4) = dt; F(2, 5) = dt; F(6, 7) = dt;
		F_arr.push_back(F);

		// Step3: Measurement Noise (R) from sigma-0
		// R = (sigma_pos)^2 * I  where sigma_pos scales with sigma-0
		// sigma0_m measures fit quality, NOT absolute accuracy.
		// Real-world SPP has systematic error the LS residual doesn't capture.
		// Inflate substantially:
		double sigma_pos = std::max(e.sigma0_m * 2.5, 15.0);
		Matrix3d R = Matrix3d::Identity() * (sigma_pos * sigma_pos);


		// Step4: Update KF with the ECEF measurement
		kf.update(z_meas, R);


		// Bias Estimator: H (3x8 position selector) -- same as inside kf.update()
		MatrixXd H = MatrixXd::Zero(3, 8);
		H(0, 0) = 1; H(1, 1) = 1; H(2, 2) = 1;

		// Bias Estimator: z_hat = H * x_pred (predicted measurement, BEFORE this update)
		Vector3d z_hat = H * x_pred_arr.back();

		// Bias Estimator: Kalman gain L used by kf.update()
		MatrixXd L = kf.getLastGain();
		if (i == 0) cout << "L dims: " << L.rows() << " x " << L.cols() << endl;
		
		biasEst.step(F, H, L, P_pred_arr.back(), R, z_meas, z_hat);
		bhat_history.push_back(biasEst.getBias());
		delta_history.push_back(biasEst.getDelta());   // <-- ADDED: save per-epoch delta
		
		// Add this debug print inside your main loop, just after delta_history.push_back(...):
		if (i == 0 || i == 100 || i == 5759) {
			cout << "  epoch " << i << "  delta = " << biasEst.getDelta().transpose() << endl;
		}

		// Step5: Extract Filtered State
		VectorXd x_est = kf.getState();
		MatrixXd P_est = kf.getCovariance();

		// SAVE filtered result (after update)
		x_filt_arr.push_back(kf.getState());
		P_filt_arr.push_back(kf.getCovariance());

		prev_t = e.t_gps;

		// Progress
		if ((i + 1) % 1000 == 0) {
			cout << "Processed " << (i + 1) << " / " << epochs.size()
				<< " epochs" << endl;
		}
	}


	// ------ RTS Smoothing -----------------------//
	//
	cout << "Running RTS backward smoother..." << endl;

	vector<VectorXd> x_smooth_arr;
	vector<MatrixXd> P_smooth_arr;

	rtsSmooth(x_filt_arr, P_filt_arr, x_pred_arr, P_pred_arr,
		F_arr, x_smooth_arr, P_smooth_arr);

	cout << "RTS smoothing complete." << endl;


	// ------ Write 4-Way Comparison CSV -----------------------//
	//
	// Variant A: kf_lat/lon/h              -- forward KF alone (Week 11 baseline)
	// Variant B: fc_lat/lon/h              -- forward KF + PER-EPOCH Friedland (rigorous)
	// Variant C: rts_lat/lon/h             -- RTS-smoothed alone (Week 12 Day 3)
	// Variant D: corrected_lat/lon/h       -- RTS-smoothed + CONVERGED-CONSTANT Friedland (approx)
	ofstream out2("../../data/trajectory_rts__crctd_Bias.csv");
	out2 << "t_gps,"
		<< "kf_lat,kf_lon,kf_h,"
		<< "fc_lat,fc_lon,fc_h,"
		<< "rts_lat,rts_lon,rts_h,"
		<< "kf_sigma,rts_sigma,"
		<< "corrected_lat,corrected_lon,corrected_h" << endl;

	for (size_t i = 0; i < epochs.size(); i++) {
		// Variant A: forward KF alone
		double kf_lat, kf_lon, kf_h;
		ecefToLla(x_filt_arr[i](0), x_filt_arr[i](1), x_filt_arr[i](2),
			kf_lat, kf_lon, kf_h);

		// Variant B (NEW): forward KF + per-epoch Friedland correction (rigorous)
		VectorXd x_filt_corr = x_filt_arr[i] + delta_history[i];
		double fc_lat, fc_lon, fc_h;
		ecefToLla(x_filt_corr(0), x_filt_corr(1), x_filt_corr(2),
			fc_lat, fc_lon, fc_h);

		// Variant C: RTS-smoothed alone
		double rts_lat, rts_lon, rts_h;
		ecefToLla(x_smooth_arr[i](0), x_smooth_arr[i](1), x_smooth_arr[i](2),
			rts_lat, rts_lon, rts_h);

		// Variant D: RTS-smoothed + converged-constant Friedland correction (approximation)
		// NOTE: biasEst.getDelta() here returns the FINAL delta (loop has already finished),
		// which is exactly what Variant D is supposed to use.
		VectorXd x_unbiased = x_smooth_arr[i] + biasEst.getDelta();
		double corrected_lat, corrected_lon, corrected_h;
		ecefToLla(x_unbiased(0), x_unbiased(1), x_unbiased(2),
			corrected_lat, corrected_lon, corrected_h);

		out2 << fixed << setprecision(6)
			<< epochs[i].t_gps << ","
			<< kf_lat << "," << kf_lon << "," << kf_h << ","
			<< fc_lat << "," << fc_lon << "," << fc_h << ","
			<< rts_lat << "," << rts_lon << "," << rts_h << ","
			<< sqrt(P_filt_arr[i](0, 0)) << "," << sqrt(P_smooth_arr[i](0, 0)) << ","
			<< corrected_lat << "," << corrected_lon << "," << corrected_h << endl;
	}
	out2.close();
	cout << "Output: trajectory_rts_crctd_Bias.csv (4-way comparison: A=kf, B=fc, C=rts, D=corrected)" << endl;

	// Print final bias estimate for a quick sanity check
	Vector3d bfinal = biasEst.getBias();
	cout << "Final bhat (ECEF, m): " << bfinal.transpose() << endl;

	return 0;
}