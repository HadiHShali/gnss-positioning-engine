// kf_spp_processor_crctd.cpp — corrected copy of kf_spp_processor.cpp
// Month 3 / Week 12 / Day 2 (review 2026-10-24)
// Pipes the single-frequency SPP trajectory (Week 9) through the Kalman filter
// and the RTS smoother -> trajectory_rts_crctd.csv
// Week 13 Day 3 (DualFreqSpp_KF_RTS_crctd.cpp) is this same file with a different
// input: `diff` the two and only the input name and labels differ.
//
// Changes vs. the original are tagged [CRCTD]. They are the Week 10 fixes, now
// with the RTS pass of Week 12 kept:
//   1. INPUT: the transmit-time-corrected Week 9 output, trajectory_crctd.csv
//      (the original read the reception-time trajectory.csv).
//   2. STATIC mode: q_accel = 0 and velocity pinned to ~0. BILL00USA does not move,
//      so letting the filter estimate a velocity turns slow SPP error drift into fake
//      motion (Week 10 finding). With the buildQ loop fixed, the old q_accel = 0.01
//      would now REALLY random-walk the velocity (~0.04 m/s per epoch).
//   3. Non-converged LS epochs are not fed into the update.
//   4. R follows each epoch's geometry and fit quality (sigma0 * PDOP) instead of a
//      constant 15 m floor that made sigma0 dead code.
//   5. Output written with 9 decimals (6 decimals of a degree = a 0.11 m x 0.09 m grid).
//   6. ADDED: 3D error summary (raw / KF / RTS vs IGS truth), the RTS offset in
//      East-North-Up, and the NEES consistency check of Week 10.
//
// STATIC-STATION NOTE: with Q = 0 the forward KF is a running weighted mean, and
// the RTS smoother gives every epoch the FINAL estimate (a flat line): the batch
// weighted mean of the whole day. Its error is the day's common offset (bias),
// which no filter can see from one station (Week 12 observability finding).

// ----------------------Section 0: Header + Setup------------------------------------------//
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>   // [CRCTD] std::max
#include "kalman_filter_gps.h"
#include "trajectory_reader.h"
#include "utils.h"

using namespace std;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::Matrix3d;


// ------ [CRCTD] Tuning constants, all in one place (same values as Week 10) ------//
// R model switch: true = per-epoch R from sigma0 * PDOP (new),
//                 false = the original constant 15 m (to compare one change at a time)
const bool   USE_GEOMETRY_R = true;
// sigma0 is 0 when n_sats == 4 (no redundancy), so it needs a floor
const double SIGMA0_FLOOR_M = 0.5;
// SPP errors (multipath, residual troposphere, orbit/clock) are correlated over
// HOURS, but the KF treats every 15 s epoch as independent. Inflating R is the
// simplest compensation. In static mode (Q = 0) a constant K changes only the
// reported sigma and NEES, NOT the position: tune it with NEES below.
const double K_INFL = 20.0;
// Skip the convergence transient when computing NEES statistics
const double NEES_START_S = 3600.0;


// [CRCTD] 3D error statistics of a set of ECEF positions against truth
struct ErrStats { double mean = 0.0, rms = 0.0, max = 0.0; };
static ErrStats errStats(const vector<Vector3d>& pos, const Vector3d& truth)
{
	ErrStats s;
	for (const auto& p : pos) {
		double e = (p - truth).norm();   // 3D distance to truth [m]
		s.mean += e;
		s.rms  += e * e;
		s.max   = std::max(s.max, e);
	}
	s.mean /= pos.size();
	s.rms   = std::sqrt(s.rms / pos.size());
	return s;
}


int main() {
	cout << "=== KF + RTS PROCESSOR, single-frequency input (crctd) ===" << endl;   // [CRCTD] title only

	// Load raw SPP trajectory
	// [CRCTD] was "../../data/trajectory.csv" (reception-time SPP).
	//   Copy trajectory_crctd.csv from the Week 9 SPP processor's data/ folder into this day's data/.
	string in_csv = "../../data/trajectory_crctd.csv";
	auto epochs = readTrajectoryCsv(in_csv);
	cout << "Loaded " << epochs.size() << " raw epochs" << endl;

	if (epochs.size() < 4) {
		cerr << "Not enough epochs!" << endl;
		return 1;
	}

	// [CRCTD] BILL00USA truth (IGS) — used only for the summary and NEES
	const double TRUE_X = -2420420.501;
	const double TRUE_Y = -4737131.551;
	const double TRUE_Z =  3507827.665;
	const Vector3d truth(TRUE_X, TRUE_Y, TRUE_Z);


	// ------ initialize the KF with the first epoch-----------------------//
	//
	// [CRCTD] static station -> NO process noise on position/velocity (was 0.01).
	KalmanFilterGps kf(0.0);

	// Initial State from FIRST epoch (should be inj ECEF)
	double x0, y0, z0;
	llaToEcef(epochs[0].lat_deg, epochs[0].lon_deg, epochs[0].height_m,
		x0, y0, z0);

	VectorXd x_init = VectorXd::Zero(8);
	x_init(0) = x0;   // position ECEF
	x_init(1) = y0;
	x_init(2) = z0;

	//Velocities and Clock terms all start at 0

	// Initialize the covariance matrix
	MatrixXd P_init = MatrixXd::Identity(8, 8);
	P_init.block(0, 0, 3, 3) *= (50.0 * 50.0);   // 50m 1-sigma pos uncertainty
	// [CRCTD] velocity is KNOWN to be ~0: sigma_v = 1e-7 m/s (< 1 cm drift over 24 h).
	//   Was 1 m/s, which let the filter fit a straight line through the error drift.
	P_init.block(3, 3, 3, 3) *= 1e-14;
	P_init(6, 6) = 100.0;                          // clock states: unobserved here (position-only measurement)
	P_init(7, 7) = 1.0;


	kf.initialize(x_init, P_init);
	cout << "KF initialized at first epoch (static mode)" << endl;


	double prev_t = epochs[0].t_gps;
	int n_skipped = 0;            // [CRCTD] non-converged epochs
	double nees_sum = 0.0;        // [CRCTD] NEES accumulator
	int nees_n = 0;


	// Storage for RTS smoother
	vector<VectorXd> x_filt_arr, x_pred_arr;
	vector<MatrixXd> P_filt_arr, P_pred_arr, F_arr;

	x_filt_arr.reserve(epochs.size());
	x_pred_arr.reserve(epochs.size());
	P_filt_arr.reserve(epochs.size());
	P_pred_arr.reserve(epochs.size());
	F_arr.reserve(epochs.size());

	vector<Vector3d> raw_ecef;    // [CRCTD] raw SPP positions (converged only), for the summary
	raw_ecef.reserve(epochs.size());

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
		if (e.converged) raw_ecef.push_back(z_meas);   // [CRCTD] summary uses real fixes only


		// Step2: Predict (skip the first epoch)
		double dt =15.0; //default for epoch 0
		if (i > 0) {
			dt = e.t_gps - prev_t;
			if (dt <= 0.0 || dt>60.0) dt = 15.0;   // safety (no effect in static mode: Q = 0)
			kf.predict(dt);
		}


		// SAVE prediction (before update)
	    x_pred_arr.push_back(kf.getState());
    	P_pred_arr.push_back(kf.getCovariance());


		// Reconstruct F used this step (matches buildF() in your class)
		MatrixXd F = MatrixXd::Identity(8, 8);
		F(0,3) = dt; F(1,4) = dt; F(2,5) = dt; F(6,7) = dt;
    	F_arr.push_back(F);

		//Step3: Measurement Noise (R)
		double sigma_pos;
		if (USE_GEOMETRY_R) {
			// [CRCTD] sigma0 * PDOP is the formal 3D position sigma of the LS fix;
			// divide by sqrt(3) for a per-axis value, then inflate for time correlation.
			double s0 = std::max(e.sigma0_m, SIGMA0_FLOOR_M);
			sigma_pos = K_INFL * s0 * e.pdop / std::sqrt(3.0);
		}
		else {
			// original model: in practice the 15 m floor won at every epoch
			sigma_pos = std::max(e.sigma0_m * 2.5, 15.0);
		}
		Matrix3d R = Matrix3d::Identity() * (sigma_pos * sigma_pos);


		// Step4: Update KF with the ECEF measurement — [CRCTD] only with converged LS fixes.
		// A failed fix carries no information: keep the prediction (filtered = predicted).
		if (e.converged) {
			kf.update(z_meas, R);
		}
		else {
			n_skipped++;
		}

		// Step5 : Extract Filtered Step
		VectorXd x_est = kf.getState();
		MatrixXd P_est = kf.getCovariance();

	    // SAVE filtered result (after update)
    	x_filt_arr.push_back(x_est);
    	P_filt_arr.push_back(P_est);

		// [CRCTD] NEES = e' P^-1 e on the position block (3 DOF -> expect ~3)
		if (e.t_gps - epochs[0].t_gps >= NEES_START_S) {
			Vector3d err = x_est.head<3>() - truth;
			Matrix3d Ppos = P_est.block<3, 3>(0, 0);
			nees_sum += err.dot(Ppos.ldlt().solve(err));
			nees_n++;
		}

		prev_t = e.t_gps;

		// Progress
		if ((i + 1) % 1000 == 0) {
			cout << "Processed " << (i + 1) << " / " << epochs.size()
				<< " epochs" << endl;
		}

	}


	// ------ RTS backward pass-----------------------//
	cout << "Running RTS backward smoother..." << endl;

	vector<VectorXd> x_smooth_arr;
	vector<MatrixXd> P_smooth_arr;

	rtsSmooth(x_filt_arr, P_filt_arr, x_pred_arr, P_pred_arr,
          	F_arr, x_smooth_arr, P_smooth_arr);

	cout << "RTS smoothing complete." << endl;

	// Write combined output CSV: KF-filtered AND RTS-smoothed
	// [CRCTD] new file name, so the original's trajectory_rts.csv is not overwritten
	ofstream out2("../../data/trajectory_rts_crctd.csv");
	out2 << "t_gps,kf_lat,kf_lon,kf_h,rts_lat,rts_lon,rts_h,"
     		<< "kf_sigma,rts_sigma" << endl;

	vector<Vector3d> kf_ecef, rts_ecef;   // [CRCTD] for the summary
	for (size_t i = 0; i < epochs.size(); i++) {
    		double kf_lat, kf_lon, kf_h;
    		ecefToLla(x_filt_arr[i](0), x_filt_arr[i](1), x_filt_arr[i](2),
              		kf_lat, kf_lon, kf_h);

    		double rts_lat, rts_lon, rts_h;
    		ecefToLla(x_smooth_arr[i](0), x_smooth_arr[i](1), x_smooth_arr[i](2),
              		rts_lat, rts_lon, rts_h);

    		kf_ecef.push_back(x_filt_arr[i].head<3>());     // [CRCTD]
    		rts_ecef.push_back(x_smooth_arr[i].head<3>());  // [CRCTD]

    		// [CRCTD] 9 decimals (was 6): 1e-9 deg ~ 0.1 mm. kf_sigma/rts_sigma = sqrt(P_xx) [m]
    		out2 << fixed << setprecision(9)
         		<< epochs[i].t_gps << ","
         		<< kf_lat << "," << kf_lon << "," << kf_h << ","
         		<< rts_lat << "," << rts_lon << "," << rts_h << ","
         		<< sqrt(P_filt_arr[i](0,0)) << ","
         		<< sqrt(P_smooth_arr[i](0,0)) << endl;
	}
	out2.close();


	// ------ [CRCTD] Summary-----------------------//
	ErrStats s_raw = errStats(raw_ecef, truth);
	ErrStats s_kf  = errStats(kf_ecef,  truth);
	ErrStats s_rts = errStats(rts_ecef, truth);

	cout << endl << "======================================" << endl;
	cout << "       KF + RTS PROCESSING COMPLETE" << endl;
	cout << "======================================" << endl;
	cout << "Total epochs processed:       " << epochs.size() << endl;
	cout << "Non-converged epochs skipped: " << n_skipped << endl;
	cout << fixed << setprecision(3);
	cout << endl << "=== 3D ERROR vs IGS TRUTH (mean / RMS / max) ===" << endl;
	cout << "Single-freq raw SPP : " << s_raw.mean << " / " << s_raw.rms << " / " << s_raw.max << " m" << endl;
	cout << "Single-freq + KF    : " << s_kf.mean  << " / " << s_kf.rms  << " / " << s_kf.max  << " m" << endl;
	cout << "Single-freq + RTS   : " << s_rts.mean << " / " << s_rts.rms << " / " << s_rts.max << " m" << endl;

	// RTS offset from truth in the local frame (at the truth position)
	double lat_t, lon_t, h_t;
	ecefToLla(TRUE_X, TRUE_Y, TRUE_Z, lat_t, lon_t, h_t);
	double phi = lat_t * M_PI / 180.0, lam = lon_t * M_PI / 180.0;
	Vector3d d = x_smooth_arr[0].head<3>() - truth;   // RTS is flat in static mode: epoch 0 = all epochs
	double dE = -sin(lam) * d(0) + cos(lam) * d(1);
	double dN = -sin(phi) * cos(lam) * d(0) - sin(phi) * sin(lam) * d(1) + cos(phi) * d(2);
	double dU =  cos(phi) * cos(lam) * d(0) + cos(phi) * sin(lam) * d(1) + sin(phi) * d(2);
	cout << "RTS offset (E, N, U): (" << dE << ", " << dN << ", " << dU << ") m"
		 << "   horizontal " << std::hypot(dE, dN) << " m" << endl;

	VectorXd xf = kf.getState();
	cout << setprecision(6)
		 << "Final velocity: (" << xf(3) << ", " << xf(4) << ", " << xf(5) << ") m/s  (should be ~0)" << endl;
	cout << setprecision(3)
		 << "Final position 1-sigma (x): " << sqrt(kf.getCovariance()(0, 0)) << " m" << endl;
	if (nees_n > 0) {
		// > 3: filter is overconfident (P too small) -> raise K_INFL
		// < 3: filter is underconfident (P too large) -> lower K_INFL
		cout << "Mean NEES after " << NEES_START_S / 3600.0 << " h: "
			 << nees_sum / nees_n << "  (expect ~3)" << endl;
	}
	cout << endl << "Output: trajectory_rts_crctd.csv" << endl;   // [CRCTD] new name

	return 0;

}
