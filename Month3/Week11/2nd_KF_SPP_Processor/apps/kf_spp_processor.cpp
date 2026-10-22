// kf_spp_processor.cpp: the executable that pipes trajectory.csv through your Kalman filter
// and then export ~5,700 smoothed positions -> trajectory_kf.csv 

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
	KalmanFilterGps kf(0.001);

	// Initial State from FIRST epoch (should be inj ECEF)
	double x0, y0, z0;
	llaToEcef(epochs[0].lat_deg, epochs[0].lon_deg, epochs[0].height_m,
		x0, y0, z0);

	VectorXd x_init = VectorXd::Zero(8);
	x_init(0) = x0;   // position ECEF
	x_init(1) = y0;
	x_init(2) = z0;

	//Velocities and Clock terms all start at 0

	// Initialize the Weight Matrix
	MatrixXd P_init = MatrixXd::Identity(8, 8);
	P_init.block(0, 0, 3, 3) *= (50.0 * 50.0);   // 50m 1-sigma pos uncertainty
	P_init.block(3, 3, 3, 3) *= (1.0 * 1.0);      // 1 m/s velocity uncertainty
	P_init(6, 6) = 100.0;                          // clock bias uncertainty
	P_init(7, 7) = 1.0;                             // clock drift uncertainty


	kf.initialize(x_init, P_init);
	cout << "KF initialized at first epoch" << endl;


	// ------ Open Output CSV file-----------------------//
	//
	ofstream out("../../data/trajectory_kf_crctd.csv");
	out << "t_gps,raw_lat,raw_lon,raw_h,kf_lat,kf_lon,kf_h,"
		<< "vx,vy,vz,sqrt_P00,sqrt_P11,sqrt_P22" << endl;

	double prev_t = epochs[0].t_gps;


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
		if (i > 0) {
			double dt = e.t_gps - prev_t;
			if (dt <= 0.0 || dt>60.0) dt = 15.0;   // safety
			kf.predict(dt);
		}


		//Step3: Measurement Noise (R) from sigma-0
		// R = (sigma_pos) ^ 2 · I  where sigma_pos scales with sigma-0
		// sigma0_m measures fit quality, NOT absolute accuracy.
		// Real-world SPP has systematic error the LS residual doesn't capture.
		// Inflate substantially — try 2-3x, or add a floor:
		double sigma_pos = std::max(e.sigma0_m * 2.5, 15.0);
		Matrix3d R = Matrix3d::Identity() * (sigma_pos * sigma_pos);


		// Step4: Update KF with the ECEF measurememtn 
		kf.update(z_meas, R);

		// Step5 : Extract Filtered Step
		VectorXd x_est = kf.getState();
		MatrixXd P_est = kf.getCovariance();


		// Convert filtered ECEF back to LLA for output
		double kf_lat, kf_lon, kf_h;
		ecefToLla(x_est(0), x_est(1), x_est(2),
			kf_lat, kf_lon, kf_h);


		// Step 6: Write CSV row
		out << fixed << setprecision(6)
			<< e.t_gps << ","
			<< e.lat_deg << "," << e.lon_deg << "," << e.height_m << ","
			<< kf_lat << "," << kf_lon << "," << kf_h << ","
			<< x_est(3) << "," << x_est(4) << "," << x_est(5) << ","
			<< sqrt(P_est(0, 0)) << "," << sqrt(P_est(1, 1)) << ","
			<< sqrt(P_est(2, 2)) << endl;


		prev_t = e.t_gps;

		// Progress
		if ((i + 1) % 1000 == 0) {
			cout << "Processed " << (i + 1) << " / " << epochs.size()
				<< " epochs" << endl;
		}

	}
	out.close();


	// ------ Summary of Output-----------------------//
	//

	cout << endl << "======================================" << endl;
	cout << "          KF PROCESSING COMPLETE" << endl;
	cout << "===========================================" << endl;
	cout << "Total epochs processed: " << epochs.size() << endl;

	VectorXd xf = kf.getState();
	MatrixXd Pf = kf.getCovariance();

	// Compare final KF position vs BILL00USA truth
	const double TRUE_X = -2420420.501;
	const double TRUE_Y = -4737131.551;
	const double TRUE_Z = 3507827.665;

	double dx = xf(0) - TRUE_X;
	double dy = xf(1) - TRUE_Y;
	double dz = xf(2) - TRUE_Z;
	double err = sqrt(dx * dx + dy * dy + dz * dz);

	cout << fixed << setprecision(3);
	cout << "Final KF ECEF: (" << xf(0) << ", " << xf(1) << ", "
		<< xf(2) << ")" << endl;
	cout << "Final 3D error: " << err << " m" << endl;
	cout << "Final velocity: (" << xf(3) << ", " << xf(4) << ", "
		<< xf(5) << ") m/s" << endl;
	cout << "Final position 1-sigma: " << sqrt(Pf(0, 0)) << " m" << endl;
	cout << endl << "Output: trajectory_kf_crctd.csv" << endl;

	return 0;

}