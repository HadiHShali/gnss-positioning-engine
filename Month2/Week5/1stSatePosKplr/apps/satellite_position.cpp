// ==========================================================// 
//						Satellite Position                   //
// =========================================================//
// Computes GPS satellite position in ECEF from broadcast ephemeris
// Implements the 9step algorithm from IS-GPS-200, Section 20.3.3.4.3

// SECTION0: Headers
#include <iostream>
#include <cmath>
#include <iomanip>
using namespace std;

// ── PHYSICAL CONSTANTS ───────────────────────────────────────────────────
const double MU_EARTH = 3.986005e14;          // Earth's GM (m^3/s^2)
const double OMEGA_EARTH = 7.2921151467e-5;      // Earth rotation (rad/s)
const double PI = 3.141592653589793;

// ── GPS EPHEMERIS DATA STRUCTURE ─────────────────────────────────────────
// All 16 broadcast parameters needed to compute satellite position
struct GpsEphemeris {
	// Reference time
	double toe;        // time of ephemeris (seconds of GPS week)

	// Six Keplerian elements
	double sqrt_a;     // square root of semi-major axis (m^0.5)
	double e;          // eccentricity (dimensionless)
	double i0;         // inclination at reference time (rad)
	double Omega0;     // right ascension at reference time (rad)
	double omega;      // argument of perigee (rad)
	double M0;         // mean anomaly at reference time (rad)

	// Perturbation rates
	double delta_n;    // mean motion correction (rad/s)
	double i_dot;      // rate of inclination (rad/s)
	double Omega_dot;  // rate of right ascension (rad/s)

	// Harmonic correction coefficients
	double Cuc, Cus;   // corrections to argument of latitude (rad)
	double Crc, Crs;   // corrections to orbit radius (m)
	double Cic, Cis;   // corrections to inclination (rad)
};

// ── SATELLITE POSITION OUTPUT ────────────────────────────────────────────
struct SatPosition {
	double X;          // ECEF X coordinate (m)
	double Y;          // ECEF Y coordinate (m)
	double Z;          // ECEF Z coordinate (m)
};

// ── HELPER: Solve Kepler's equation iteratively ──────────────────────────
// Given mean anomaly M and eccentricity e, find eccentric anomaly E
// such that:  M = E - e * sin(E)
double solveKepler(double M, double e) {
	double E = M;  // initial guess
	const double TOLERANCE = 1e-12;
	const int MAX_ITER = 20;

	for (int i = 0; i < MAX_ITER; i++) {
		double E_new = M + e * sin(E);
		if (abs(E_new - E) < TOLERANCE) {
			return E_new;  // converged!
		}
		E = E_new;
	}

	// Convergence not achieved (extremely rare for GPS orbits)
	cerr << "Warning: Kepler iteration did not converge." << endl;
	return E;
}


// MAin Computation - 9step algorithm
// This function gives the X Y Z coordinates of the Satellite posotion.
// Return by SatPosition struct — output multiple values (X, Y, Z) in one bundle
// ── MAIN COMPUTATION: 9-step algorithm ───────────────────────────────────
SatPosition computeSatPosECEF(const GpsEphemeris& eph, double t_gps) {
	SatPosition pos;

	// Step 1: Compute mean motion n
	double a = eph.sqrt_a * eph.sqrt_a;       // semi-major axis
	double n0 = sqrt(MU_EARTH / (a * a * a));     // nominal mean motion
	double n = n0 + eph.delta_n;              // corrected mean motion

	// Step 2: Time since reference epoch
	double tk = t_gps - eph.toe;
	// Handle GPS week rollover
	if (tk > 302400.0) tk -= 604800.0;
	if (tk < -302400.0) tk += 604800.0;

	// Step 3: Mean anomaly at time t
	double Mk = eph.M0 + n * tk;

	// Step 4: Solve Kepler's equation for eccentric anomaly
	double Ek = solveKepler(Mk, eph.e);

	// Step 5: True anomaly
	double sin_E = sin(Ek);
	double cos_E = cos(Ek);
	double nuk = atan2(sqrt(1 - eph.e * eph.e) * sin_E, cos_E - eph.e);

	// Step 6: Argument of latitude with corrections
	double phi_k = nuk + eph.omega;
	double sin_2phi = sin(2 * phi_k);
	double cos_2phi = cos(2 * phi_k);

	double u_k = phi_k + eph.Cuc * cos_2phi + eph.Cus * sin_2phi;
	double r_k = a * (1 - eph.e * cos_E) + eph.Crc * cos_2phi + eph.Crs * sin_2phi;
	double i_k = eph.i0 + eph.i_dot * tk + eph.Cic * cos_2phi + eph.Cis * sin_2phi;

	// Step 7: Position in orbital plane
	double x_orb = r_k * cos(u_k);
	double y_orb = r_k * sin(u_k);

	// Step 8: Corrected longitude of ascending node
	double Omega_k = eph.Omega0 + (eph.Omega_dot - OMEGA_EARTH) * tk
		- OMEGA_EARTH * eph.toe;

	// Step 9: Rotate to ECEF
	double cos_Omega = cos(Omega_k);
	double sin_Omega = sin(Omega_k);
	double cos_i = cos(i_k);
	double sin_i = sin(i_k);

	pos.X = x_orb * cos_Omega - y_orb * cos_i * sin_Omega;
	pos.Y = x_orb * sin_Omega + y_orb * cos_i * cos_Omega;
	pos.Z = y_orb * sin_i;

	return pos;
}

// ── MAIN — Test with a real GPS ephemeris ────────────────────────────────
int main() {
	cout << "=== GPS Satellite Position Calculator ===" << endl;
	cout << "Implementing IS-GPS-200 Section 20.3.3.4.3" << endl << endl;

	// ── REAL EPHEMERIS DATA ─────────────────────────────────────────────
	// GPS PRN 5, January 7 2024, ~00:00 UTC
	// (Pulled from a real BRDC navigation file)
	GpsEphemeris eph;
	eph.toe = 7200.0;                    // 02:00:00 GPS time
	eph.sqrt_a = 5153.65040588379;          // ~26,560 km orbit
	eph.e = 0.00578925572569;
	eph.i0 = 0.96218657353945;          // ~55°
	eph.Omega0 = -2.46893573783000;
	eph.omega = 0.93755649238080;
	eph.M0 = -2.07477835908625;
	eph.delta_n = 4.42456802576e-9;
	eph.i_dot = -2.0608988e-10;
	eph.Omega_dot = -8.0033236687e-9;
	eph.Cuc = -1.42537057399e-6;
	eph.Cus = 9.50321555138e-6;
	eph.Crc = 199.40625000000;
	eph.Crs = -26.50000000000;
	eph.Cic = -1.490116119e-8;
	eph.Cis = -1.564621925e-7;

	// Compute position at ephemeris reference time (tk = 0)
	double t_gps = eph.toe;
	SatPosition pos = computeSatPosECEF(eph, t_gps);

	cout << fixed << setprecision(4);
	cout << "PRN 5 position at toe (" << t_gps << " s of week):" << endl;
	cout << "  X = " << pos.X << " m" << endl;
	cout << "  Y = " << pos.Y << " m" << endl;
	cout << "  Z = " << pos.Z << " m" << endl;

	// Sanity check: magnitude should be ~26,560 km (GPS orbit radius)
	double r = sqrt(pos.X * pos.X + pos.Y * pos.Y + pos.Z * pos.Z);
	cout << endl << "Distance from Earth center: " << r / 1000.0 << " km" << endl;
	cout << "Expected: ~26,560 km (GPS orbit radius)" << endl;

	return 0;
}
