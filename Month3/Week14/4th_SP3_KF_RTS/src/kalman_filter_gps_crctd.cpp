// =====================================================================//
// kalman_filter_gps_crctd.cpp — corrected copy of kalman_filter_gps.cpp //
// =====================================================================//
// Month 3 / Week 13 / Day 3 (review 2026-10-24). The class part is IDENTICAL to
// the Week 10 kalman_filter_gps_crctd.cpp; rtsSmooth() (added in Week 12) is at the end.
// Changes vs. the original are tagged [CRCTD]:
//   1. buildQ(): loop counter `i` was never initialized -> undefined behavior.
//      In your build the loop never ran, so Q had ZERO position/velocity noise
//      no matter what q_accel you passed. That is why the KF sigma curve
//      matched a Q = 0 model exactly.
//   2. update(): Joseph-form covariance update (stays symmetric and positive
//      definite over thousands of updates; the short form (I-KH)P does not).
//   3. Constructor initializer list reordered to match the declaration order
//      in the header (silences -Wreorder; no behavior change).
//   4. Added <stdexcept> (std::runtime_error was compiling only by luck,
//      through Eigen's own includes).
//   rtsSmooth() itself is unchanged (a static-mode note was added above it).

#include "kalman_filter_gps.h"
#include <stdexcept>   // [CRCTD] std::runtime_error lives here

using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::Vector3d;
using Eigen::Matrix3d;


// [CRCTD] Members are initialized in the order they are DECLARED in the header
// (x_, P_, ..., q_accel_, q_clk_bias_, q_clk_drift_), so the list now follows that order.
KalmanFilterGps::KalmanFilterGps(double q_accel,
	double q_clk_bias,
	double q_clk_drift)
	: x_(VectorXd::Zero(8)),                 // state starts at zero
	P_(MatrixXd::Identity(8, 8) * 1e6),      // huge initial uncertainty until initialize()
	q_accel_(q_accel),                       // see the units note in buildQ()
	q_clk_bias_(q_clk_bias),
	q_clk_drift_(q_clk_drift) {
}


// The '::' means: the function 'initialize' belongs to the class 'KalmanFilterGps'
void KalmanFilterGps::initialize(const VectorXd& x0, const MatrixXd& P0)
{
	if (x0.size() != 8 || P0.rows() != 8 || P0.cols() != 8)
	{
		throw std::runtime_error("KF initialize: dimensions must be 8");
	}
	x_ = x0;            // copy the starting state
	P_ = P0;            // copy the starting covariance
	initialized_ = true;
}

// --------- buildF() helper ----------------//
MatrixXd KalmanFilterGps::buildF(double dt) const
{
	// 8x8 identity + dt in position-velocity + clock-drift couplings
	MatrixXd F = MatrixXd::Identity(8, 8);
	F(0, 3) = dt;  // x += vx * dt
	F(1, 4) = dt;  // y += vy * dt
	F(2, 5) = dt;  // z += vz * dt
	F(6, 7) = dt;  // dt_bias += dt_drift * dt
	return F;
}

// --------- buildQ() helper ----------------//
// Continuous white-noise-acceleration model, per axis:
//     Q_axis = q * [ dt^3/3   dt^2/2 ]
//                  [ dt^2/2   dt     ]
// where q is a power spectral density in m^2/s^3.
// UNITS NOTE: the code uses q = q_accel^2, so q_accel is sqrt(PSD) in m/s^1.5,
// NOT an acceleration standard deviation in m/s^2. Example: q_accel = 0.001
// lets velocity random-walk by sqrt(1e-6 * 15) ~ 0.004 m/s per 15 s epoch.
// For a STATIC station the correct value is q_accel = 0.
MatrixXd KalmanFilterGps::buildQ(double dt) const
{
	double dt2 = dt * dt;               // dt^2
	double dt3 = dt2 * dt;              // dt^3
	double q_a_sq = q_accel_ * q_accel_; // PSD q (see units note above)

	// 8x8 zero matrix
	MatrixXd Q = MatrixXd::Zero(8, 8);

	// position and velocity blocks (3 axes: x, y, z)
	// [CRCTD] was `for (int i; ...)` — i held stack garbage, so the loop was skipped
	for (int i = 0; i < 3; i++)
	{
		Q(i, i) = q_a_sq * dt3 / 3.0;          // pos-pos
		Q(i + 3, i + 3) = q_a_sq * dt;         // vel-vel
		Q(i, i + 3) = q_a_sq * dt2 / 2.0;      // pos-vel cross
		Q(i + 3, i) = q_a_sq * dt2 / 2.0;      // vel-pos cross
	}
	// clock bias and drift block (states 6 and 7)
	// NOTE: in the SPP-position filter (kf_spp_processor) nothing observes these
	// two states, so their variance just grows. Harmless (no coupling to
	// position), but they carry no information there.
	Q(6, 6) = q_clk_bias_ * dt + q_clk_drift_ * dt3 / 3.0;
	Q(7, 7) = q_clk_drift_ * dt;
	Q(6, 7) = q_clk_drift_ * dt2 / 2.0;
	Q(7, 6) = q_clk_drift_ * dt2 / 2.0;
	return Q;
}

// --------- predict() ----------------//
void KalmanFilterGps::predict(double dt)
{
	if (!initialized_)
	{
		throw std::runtime_error("KF predict: not initialized");
	}

	MatrixXd F = buildF(dt);   // state transition
	MatrixXd Q = buildQ(dt);   // process noise

	// the two prediction equations:
	x_ = F * x_;                          // propagate the state
	P_ = F * P_ * F.transpose() + Q;      // propagate the uncertainty
}

// --------- update() ----------------//
void KalmanFilterGps::update(const Vector3d& z, const Matrix3d& R)
{
	if (!initialized_)
	{
		throw std::runtime_error("KF update: not initialized");
	}

	// H = Observation Matrix (3x8)
	// Measurement is z = [x, y, z] — first 3 states directly observed
	MatrixXd H = MatrixXd::Zero(3, 8);
	H(0, 0) = 1.0; // z_x = state_x
	H(1, 1) = 1.0; // z_y = state_y
	H(2, 2) = 1.0; // z_z = state_z

	// Innovation (measurement residual)
	Vector3d y = z - H * x_;

	// Innovation covariance
	Matrix3d S = H * P_ * H.transpose() + R;

	// Kalman gain (solve instead of an explicit inverse: same result, better conditioned)
	MatrixXd K = (S.ldlt().solve(H * P_)).transpose();   // = P H' S^-1 because P and S are symmetric

	// state update
	x_ = x_ + K * y;

	// [CRCTD] Joseph-form covariance update:
	//   P = (I - KH) P (I - KH)' + K R K'
	// Algebraically equal to (I - KH) P, but it cannot lose symmetry or
	// positive-definiteness through round-off.
	MatrixXd I_KH = MatrixXd::Identity(8, 8) - K * H;
	P_ = I_KH * P_ * I_KH.transpose() + K * R * K.transpose();
	P_ = 0.5 * (P_ + P_.transpose());   // force exact symmetry
}

// -------------------RTS Smoother------------------------//
// Rauch-Tung-Striebel: after the forward KF pass, sweep BACKWARD so every epoch
// is estimated from ALL the data (past and future).
//
// STATIC-STATION NOTE (Q = 0, velocity pinned to ~0):
//   then P_pred[k+1] = F P_filt[k] F', so the gain C_k = P_filt F' P_pred^-1 = F^-1,
//   and x_smooth[k] = F^-1 x_smooth[k+1]. With v ~ 0, F^-1 changes nothing:
//   EVERY epoch gets the final, all-data estimate. For a fixed antenna the RTS
//   result is the batch weighted mean of the whole day, a flat line.
void rtsSmooth(
    const std::vector<VectorXd>& x_filt,
    const std::vector<MatrixXd>& P_filt,
    const std::vector<VectorXd>& x_pred,
    const std::vector<MatrixXd>& P_pred,
    const std::vector<MatrixXd>& F_used,
    std::vector<VectorXd>& x_smooth,
    std::vector<MatrixXd>& P_smooth)
{
    int N = static_cast<int>(x_filt.size());
    x_smooth.resize(N);
    P_smooth.resize(N);

    // Last epoch: smoothed = filtered (no future data available)
    x_smooth[N - 1] = x_filt[N - 1];
    P_smooth[N - 1] = P_filt[N - 1];

    // Sweep backward: k = N-2 down to 0
    for (int k = N - 2; k >= 0; --k) {
        // F_used[k+1] is the F that was used to go FROM epoch k
        //   TO the prediction x_pred[k+1]
        const MatrixXd& F = F_used[k + 1];

        // Smoother gain: C_k = P_filt[k] * F^T * inv(P_pred[k+1])
        MatrixXd Ck = P_filt[k] * F.transpose() * P_pred[k + 1].inverse();

        // Smoothed state
        x_smooth[k] = x_filt[k] +
            Ck * (x_smooth[k + 1] - x_pred[k + 1]);

        // Smoothed covariance
        P_smooth[k] = P_filt[k] +
            Ck * (P_smooth[k + 1] - P_pred[k + 1]) * Ck.transpose();
    }
}
