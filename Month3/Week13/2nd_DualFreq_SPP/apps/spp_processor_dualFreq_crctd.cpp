// spp_processor_dualFreq_crctd.cpp - Dual-frequency (ionosphere-free) SPP: every epoch -> trajectory CSV
//
// CORRECTED copy of spp_processor_dualFreq.cpp (Month 3 / Week 13 review, 2026-10-24).
// The original file stays untouched, so `diff` between the two shows exactly what
// changed. Every change is marked with "CORRECTED".
//
// WHAT THE ORIGINAL ALREADY HAD RIGHT (dual-frequency part, the point of Week 13):
//   - ionosphere-free code P_IF = K1*P1 - K2*P2 from C1W/C2W, with Klobuchar skipped for
//     those satellites (the combination already removes ~99.9% of the ionosphere)
//   - TGD skipped for P_IF and applied for single-frequency P1. IS-GPS-200: the broadcast
//     clock (af0, af1, af2) is estimated by the Control Segment from the L1/L2 P(Y)
//     ionosphere-free combination, so TGD cancels in P_IF and is needed only for L1 alone.
//   - satellite clock with af2, the relativistic term and the week wrap; troposphere for all
//   - 3D error against BILL's IGS position
//
// WHAT WAS STILL WRONG - the same timing bug as Weeks 6-9, inside the epoch loop:
//   The receiver tags each epoch with the RECEPTION time t_obs (epoch.t_gps). The signal
//   left the satellite ~67-86 ms earlier, at the TRANSMIT time t_tx. In the Earth-fixed
//   frame a GPS satellite moves ~3 km/s, so in that time it travels ~200-300 m.
//     1. Satellite position : computeSatPosECEF(eph, t_obs)            -> now at t_tx
//     2. Earth rotation     : travel time = P_IF/c (contains BOTH clocks) -> now the true
//                             travel time (t_obs - dt_rx) - t_tx
//     3. Satellite clock    : evaluated at t_obs, referenced to toe      -> now at t_tx, from toc
//   Dual frequency does NOT hide this bug: P_IF removes the ionosphere, not geometry.
//
// SMALLER CORRECTIONS (no effect on the position):
//     4. n_iono_free / n_single_freq were counted BEFORE the 10-degree mask, so they did
//        not add up to n_sats (counted after it). Now counted on the satellites actually used.
//     5. CSV lat/lon had 6 decimals = a 0.11 m (north) x 0.09 m (east) grid, the staircase
//        seen in the Week 10 KF plots. Now 9 decimals (0.1 mm). Same fix as Week 10.
//     6. The debug print set cout precision to 3 and only unset 'fixed', so precision 3
//        stayed on: "Success rate" printed 3 significant digits (96.28% -> "96.3%",
//        99.98% -> "100%"). The old precision is now restored.
//     7. Output file name: trajectory_dualfreq_crctd.csv (original output not overwritten).
//
// HOW THE FIX FITS EACH EPOCH (identical to Week 9 _crctd):
//   satellite state at t_tx, ONCE per satellite per epoch (see computeSatGeom)
//   -> build obs with dt_rx = 0 -> rough solve (position AND receiver clock dt_rx)
//   -> REBUILD obs with the rough dt_rx (Sagnac now uses the true flight time)
//   -> mask -> atmosphere -> final solve   (unchanged from here on)
//
// NOT CHANGED: satellites.cpp, rinex_obs.cpp, gps_ls_solver.cpp, utils.cpp.
//   The helper functions below are copied unchanged from Week 9 _crctd
//   (= Week 6 Day 4 = Week 7 = Week 8 = Week 16 Day 3 sat_geometry.cpp math).
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "gps_ls_solver.h"
#include "utils.h"
#include "satellites.h"
#include "rinex_obs.h"

using namespace std;
using namespace Eigen;

const double C_LIGHT = 299792458.0;
const double OMEGA_EARTH = 7.2921151467e-5;

// ── CORRECTED: transmit-time helpers (copied unchanged from Week 9 _crctd) ──
// Everything inside "namespace { }" is private to this file, like 'static'.
namespace {

const double MU_GPS    = 3.986005e14;       // Earth GM used by GPS (IS-GPS-200) [m^3/s^2]; same value as satellites.cpp
const double F_REL     = -4.442807633e-10;  // relativistic clock constant -2*sqrt(MU)/c^2 [s/sqrt(m)] (IS-GPS-200)
const double HALF_WEEK = 302400.0;          // half a GPS week [s]
const double WEEK      = 604800.0;          // one GPS week [s]

// Wrap a time DIFFERENCE into [-half week, +half week].
//
// WHY: all GPS times here are "seconds of week" (0 ... 604800 s). The counter
// RESETS to 0 every Sunday 00:00 GPS time. A difference that crosses the reset
// comes out wrong by one full week, unless we fix it.
//
// ANALOGY - a 12-hour clock: it is 1 o'clock, something happened at 11 o'clock.
//   Naive: 1 - 11 = -10 h.  Real: 2 h ago.  Fix: add one full cycle: -10 + 12 = +2 h.
// Here the "cycle" is one GPS week (604800 s):
//
//   Case                      t                    toc                  naive t-toc   wrapped    meaning
//   just AFTER Sunday 00:00   0                    597600 (Sat 22:00,   -597600 s     +7200 s    toc was 2 h ago
//                                                   last week)
//   just BEFORE Sunday 00:00  604000 (Sat 23:46)   0 (next Sun 00:00)   +604000 s     -800 s     toc is 13 min ahead
//
// WHY THE THRESHOLD IS HALF A WEEK (302400 s): a valid ephemeris is only used a few
//   hours from its toe/toc (fit interval ~4 h). A real difference is never anywhere
//   near 3.5 days, so anything bigger can only mean the counter wrapped -> add or
//   subtract one week.
//
// WHY IT MATTERS: dt multiplies the clock drift af1. Without the wrap, a typical
//   af1 ~ 1e-11 s/s times 597600 s ~ 6 microseconds ~ 1.8 km of range error.
//   (In eccentricAnomaly, tk would also propagate the orbit a whole week off.)
//
// SAME LOGIC ALREADY EXISTS in satellites.cpp, computeSatPosECEF(), lines 151-152,
//   which wraps tk = t - toe the same way (and in the original's clock block, lines 106-107).
//
// WHEN IT TRIGGERS: only for epochs near Saturday/Sunday midnight. This file processes
//   the whole of Sunday 2024-01-07, so the first epochs after 00:00 GPS can use ephemerides
//   whose toc/toe are still Saturday's: there the wrap DOES matter.
double wrapWeek(double dt)
{
    if (dt >  HALF_WEEK) dt -= WEEK;   // t at week end, toc in next week:  604000 - 0      -> -800 s
    if (dt < -HALF_WEEK) dt += WEEK;   // t at week start, toc in last week:     0 - 597600 -> +7200 s
    return dt;
}

// Eccentric anomaly E_k at time t [rad]. Needed only for the relativistic clock term.
// Same first steps as computeSatPosECEF(): mean motion, mean anomaly, Kepler's equation.
double eccentricAnomaly(const GpsEphemeris& eph, double t)
{
    double a  = eph.sqrt_a * eph.sqrt_a;                     // semi-major axis [m]
    double n  = sqrt(MU_GPS / (a * a * a)) + eph.delta_n;    // corrected mean motion [rad/s]
    double tk = wrapWeek(t - eph.toe);                       // time since toe [s]
    double M  = eph.M0 + n * tk;                             // mean anomaly [rad]
    double E  = M;                                           // start value (e is small, so E ~ M)
    for (int i = 0; i < 10; ++i)                             // Newton's method on E - e*sin(E) - M = 0
        E -= (E - eph.e * sin(E) - M) / (1.0 - eph.e * cos(E));
    return E;
}

// Days since 1970-01-01 for a calendar date (H. Hinnant's days_from_civil).
long daysFromCivil(int y, int m, int d)
{
    y -= (m <= 2);                                           // Jan/Feb count as months 13/14 of the previous year
    const long era = (y >= 0 ? y : y - 399) / 400;           // 400-year cycle number
    const long yoe = y - era * 400;                          // year of era [0, 399]
    const long doy = (153 * (m > 2 ? m - 3 : m + 9) + 2) / 5 + d - 1;   // day of year, counted from March 1
    const long doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;  // day of era
    return era * 146097 + doe - 719468;                      // shift so that 1970-01-01 = day 0
}

// Time of clock (toc) as GPS seconds of week [s].
// toc = "time of clock": the reference epoch of the clock polynomial
//   dt_sat = af0 + af1*(t - toc) + af2*(t - toc)^2
// (toe = "time of ephemeris" is the reference epoch of the ORBIT parameters instead.)
// toc is the date/time on the FIRST line of the nav record (eph.year ... eph.second),
// e.g. "G05 2024 01 07 12 00 00", so it must be converted to seconds of week here.
// IS-GPS-200 references the clock polynomial (af0, af1, af2) to toc, NOT to toe.
// For GPS the two are normally equal, but toc is the correct reference.
double tocSecondsOfWeek(const GpsEphemeris& eph)
{
    int year = eph.year;
    if (year < 100) year += (year < 80) ? 2000 : 1900;       // RINEX 2 writes 2-digit years

    // GPS time started on Sunday 1980-01-06 00:00. Days since then, modulo 7,
    // is the day of the GPS week (0 = Sunday).
    long days = daysFromCivil(year, eph.month, eph.day) - daysFromCivil(1980, 1, 6);
    long dow  = ((days % 7) + 7) % 7;                        // keep it in 0..6

    return dow * 86400.0 + eph.hour * 3600.0 + eph.minute * 60.0 + eph.second;
}

// Everything this app needs about ONE satellite at ONE epoch.
struct SatGeom {
    double t_tx   = 0.0;   // signal transmit time, GPS seconds of week [s]
    double dt_sat = 0.0;   // satellite clock offset at t_tx [s] = polynomial + relativistic (TGD applied in main)
    double dt_rel = 0.0;   // relativistic part of dt_sat [s] (kept separately for debugging prints)
    double X = 0.0;        // satellite ECEF position at t_tx [m],
    double Y = 0.0;        //   still in the Earth-fixed frame OF t_tx
    double Z = 0.0;        //   (Earth rotation during the flight is applied in main)
};

// Transmit time, satellite clock and satellite position from ONE pseudorange.
//   eph   ephemeris for this satellite (from findBestEphemeris)
//   t_obs receiver time tag of the epoch [s of week] (epoch.t_gps)
//   P     measured pseudorange [m]
SatGeom computeSatGeom(const GpsEphemeris& eph, double t_obs, double P)
{
    SatGeom g;

    // 1. Pseudorange-based transmit time.
    //    P     = c*(t_rx_true - t_tx_true) + c*dt_rx - c*dt_sat    (+ atmosphere, noise)
    //    t_obs = t_rx_true + dt_rx                                  (receiver clock reading)
    //    =>  t_obs - P/c = t_tx_true + dt_sat     -> the receiver clock dt_rx CANCELS,
    //    so the transmit time is known BEFORE solving for the receiver.
    //    Atmosphere (~10 m = 33 ns) moves the satellite only ~0.1 mm: negligible.
    const double t_raw = t_obs - P / C_LIGHT;

    // 2. Satellite clock at the transmit time. dt_sat needs t_tx, and t_tx needs dt_sat,
    //    so iterate twice (the second pass changes the result by < 1e-15 s).
    const double toc = tocSecondsOfWeek(eph);
    double dt_sat = 0.0, dt_rel = 0.0;
    for (int k = 0; k < 2; ++k) {
        double t    = t_raw - dt_sat;                        // current transmit-time estimate [s]
        double dtc  = wrapWeek(t - toc);                     // time since toc [s], week-wrapped
        double poly = eph.clk_bias                           // af0 [s]
                    + eph.clk_drift      * dtc               // af1 [s/s]
                    + eph.clk_drift_rate * dtc * dtc;        // af2 [s/s^2]
        dt_rel = F_REL * eph.e * eph.sqrt_a                  // relativistic term (orbit eccentricity),
               * sin(eccentricAnomaly(eph, t));              //   amplitude ~23 ns (~7 m) per 0.01 of e
        dt_sat = poly + dt_rel;
    }
    g.dt_sat = dt_sat;
    g.dt_rel = dt_rel;

    // 3. True transmit time, and the satellite position AT THAT TIME.
    g.t_tx = t_raw - dt_sat;
    SatPosition p = computeSatPosECEF(eph, g.t_tx);
    g.X = p.X;  g.Y = p.Y;  g.Z = p.Z;
    return g;
}

// CORRECTED: one satellite's inputs for one epoch, computed once and reused by both solves.
// Same as Week 9 _crctd plus 'iono_free', which tells buildObservations whether TGD applies.
struct SatInput {
    int prn;          // satellite PRN
    double P;         // pseudorange used in the solve [m]: P_IF (dual) or P1 = C1W (single)
    bool iono_free;   // true -> P is P_IF: no TGD, no Klobuchar
    SatGeom g;        // transmit time, satellite clock, satellite position at t_tx
    double tgd;       // L1/L2 group delay from the nav record [s]
};

// CORRECTED: build the solver's Observation vector for one epoch, for a given receiver clock.
// The Earth-rotation (Sagnac) correction needs the TRUE travel time:
//   travel time = t_rx_true - t_tx = (t_obs - dt_rx) - t_tx
// The original used P_IF/c, which also contains c*dt_rx and -c*dt_sat. Each millisecond
// of clock error rotates the satellite ~2 m too far (OMEGA_EARTH * 1e-3 s * 26,560 km),
// which moved the position ~0.4 m per ms in the Week 6 test.
// dt_rx is unknown until we solve, so this is called twice per epoch: dt_rx = 0 for the
// rough solve, then the rough solve's clock estimate for everything after it.
vector<Observation> buildObservations(const vector<SatInput>& sats, double t_obs, double dt_rx)
{
    vector<Observation> out;
    for (const auto& s : sats) {
        // Sagnac
        // CORRECTED: was  travel_time = base_pr / C_LIGHT
        double travel_time = (t_obs - dt_rx) - s.g.t_tx;            // true signal flight time [s]
        double rot = OMEGA_EARTH * travel_time;                     // Earth's rotation during the flight [rad]
        double cr = cos(rot), sr = sin(rot);
        // Rotate the satellite from the Earth-fixed frame of t_tx into the frame of
        // the reception time (same rotation as the original, only the angle changed).
        double sx =  cr * s.g.X + sr * s.g.Y;
        double sy = -sr * s.g.X + cr * s.g.Y;

        Observation o;
        o.prn = s.prn;
        o.sat_x = sx; o.sat_y = sy; o.sat_z = s.g.Z;                // rotation about Z leaves Z unchanged
        // Model: P + c*dt_sat(signal) = rho + c*dt_rx  (what the solver expects)
        // CORRECTED: dt_sat now at t_tx, from toc.
        // TGD as in the original (IS-GPS-200 20.3.3.3.3.2-3):
        //   P_IF : dt_sat(IF) = dt_sat           (broadcast clock IS the P1/P2 iono-free clock)
        //   P1   : dt_sat(L1) = dt_sat - TGD
        double dt_sig = s.iono_free ? s.g.dt_sat : (s.g.dt_sat - s.tgd);
        o.pseudorange = s.P + C_LIGHT * dt_sig;
        o.used_iono_free = s.iono_free;                             // Step 4 skips Klobuchar for these
        out.push_back(o);
    }
    return out;
}

} // namespace
// ─────────────────────────────────────────────────────────────────────────────

int main()
{
	cout << "=== FULL SPP PROCESSOR (dual-frequency, transmit-time corrected) ===" << endl;   // CORRECTED: title only

	// ---------------------------------------------------------------------------------//
	//								One_time Step (before loop)							//
	// ---------------------------------------------------------------------------------//
    string nav_file = "../../data/BRDC00IGS_R_20240070000_01D_MN.rnx";
    string obs_file = "../../data/BILL00USA_R_20240070000_01D_15S_MO.rnx";

	// Parse Nav file
    vector<GpsEphemeris> ephs = parseRinexNav(nav_file);
	KlobucharParams kp = parseKlobucharFromNav(nav_file);
    cout << "Parsed " << ephs.size() << " ephemerides" << endl;

	// Parse All observation epochs
	vector<ObsEpoch> epochs = parseRinexAllEpochs(obs_file);
	cout << "Parsed " << epochs.size() << " epochs" << endl;

	// Open trajectoory CSV
	// CORRECTED (output file name only): was "../../data/trajectory_dualfreq.csv".
	//   A separate name keeps the original's trajectory from being overwritten.
	ofstream csv("../../data/trajectory_dualfreq_crctd.csv");
	csv << "t_gps,lat_deg,lon_deg,height_m,sigma0_m,pdop,n_sats,converged,n_iono_free,n_single_freq,err_3d" << endl;

	// Track statistics
	int successful = 0;
	int failed = 0;
	// NEW: for 3D error summary
	const double TRUE_X = -2420420.5012;
	const double TRUE_Y = -4737131.5508;
	const double TRUE_Z =  3507827.665;
	vector<double> errors_3d;   // store every successful epoch's error
	// ---------------------------------------------------------------------------------//
	//								Loop. Process every Epoch							//
	// ---------------------------------------------------------------------------------//

	for (const auto& epoch: epochs)
	{
		// skip epochs with too few satellites
		if (epoch.gps_records.size() < 4){failed++; continue;}

		// Step1: CORRECTED — satellite state at the TRANSMIT time, once per satellite
		// (sat clock + relativity at t_tx from toc; TGD only for single frequency;
		//  Sagnac applied in buildObservations)
		vector<SatInput> sats;
		for (const auto& rec: epoch.gps_records)
		{
			// Picking the ephemeris with t_obs is fine: 75 ms never changes which toe is nearest.
			const GpsEphemeris* eph = findBestEphemeris(ephs, rec.prn, epoch.t_gps);
			if (!eph) continue;

			// Choose pseudorange source based on C2W availability (unchanged)
			double base_pr;
			bool is_iono_free=false;
			if (rec.has_c2){
				base_pr = ionoFreePseudorange(rec.pseudorange, rec.pseudorange_c2);
				is_iono_free = true;
			} else{
				base_pr = rec.pseudorange;
				is_iono_free = false;
			}
			// CORRECTED: the per-epoch counters moved to Step 4 (after the elevation mask),
			//            so n_iono_free + n_single_freq = n_sats in the CSV.

			    // ────────── ADD DEBUG PRINT RIGHT HERE ──────────
			static bool printed_once = false;
			if (!printed_once && is_iono_free) {
				// CORRECTED: remember the current precision, so it can be put back below.
				//   unsetf(ios::fixed) resets only the 'fixed' flag; precision 3 stayed on and
				//   the "Success rate" line at the end printed 3 significant digits.
				streamsize old_prec = cout.precision();
				cout << fixed << setprecision(3);   // ADD THIS — show 3 decimal places
				cout << "  DEBUG sat G" << rec.prn
					 << "  P1(C1W)=" << rec.pseudorange
					 << "  P2(C2W)=" << rec.pseudorange_c2
					 << "  P_IF="    << base_pr
					 << "  (P2-P1)=" << (rec.pseudorange_c2 - rec.pseudorange)
					 << endl;
				cout.unsetf(ios::fixed);   // reset formatting for the rest of the program
				cout.precision(old_prec);  // CORRECTED: also reset the precision
				printed_once = true;
			}
			// ──────────────────────────────────────────────

			// CORRECTED: was  SatPosition sp = computeSatPosECEF(*eph, epoch.t_gps);   (RECEPTION time)
			// CORRECTED: the original's clock block (af0 + af1*dt + af2*dt^2 + relativistic,
			//            week-wrapped) was evaluated at the RECEPTION time and referenced to toe.
			//            computeSatGeom() evaluates the same terms at the TRANSMIT time, from toc.
			// The transmit time comes from the measured P1 (C1W). P_IF would give the same
			// t_tx to ~0.1 mm of satellite motion (the ionosphere is ~33 ns), but has ~3x the noise.
			SatInput s;
			s.prn       = rec.prn;
			s.P         = base_pr;                 // P_IF or P1: what the solve uses
			s.iono_free = is_iono_free;
			s.g         = computeSatGeom(*eph, epoch.t_gps, rec.pseudorange);
			s.tgd       = eph->TGD;                // applied only when !iono_free (buildObservations)
			sats.push_back(s);
		}

		// CORRECTED: the rough solve does not know the receiver clock yet, so it uses dt_rx = 0.
		vector<Observation> obs = buildObservations(sats, epoch.t_gps, 0.0);
		if (obs.size() < 4) { failed++; continue; }


		// Step2: Rough Solve
		        // ── STEP 2: Rough solve ──
		// (rough position for the mask and the atmosphere, and, CORRECTED, the receiver clock)
        PositionSolution rough = solveGpsPosition(obs, false);
        if (!rough.converged) { failed++; continue; }

		// CORRECTED: rebuild with the rough solve's receiver clock, so the Sagnac rotation
		// uses the true flight time. Satellite times and clocks don't change; only the
		// rotation angle does (by OMEGA_EARTH * dt_rx). Everything below uses these.
		obs = buildObservations(sats, epoch.t_gps, rough.clock_bias);


// ── STEP 3: Elevation mask ──
        vector<Observation> masked;
        for (const auto& o : obs) {
            double el = computeElevation(rough.x, rough.y, rough.z,
                                          o.sat_x, o.sat_y, o.sat_z);
            if (el >= 10.0) masked.push_back(o);
        }
        if (masked.size() < 4) { failed++; continue; }



		double rough_lat, rough_lon, rough_h;
        ecefToLla(rough.x, rough.y, rough.z, rough_lat, rough_lon, rough_h);


        // ── STEP 4: Apply atmospheric corrections ──
        // (unchanged: subtracted, P - I - T; Klobuchar only for single-frequency satellites)
		int n_iono_free_epoch = 0, n_single_freq_epoch = 0;   // CORRECTED: counted here, after the mask
        vector<Observation> corrected;
        for (const auto& o : masked) {
            double el = computeElevation(rough.x, rough.y, rough.z,
                                          o.sat_x, o.sat_y, o.sat_z);
            double az = computeAzimuth(rough.x, rough.y, rough.z,
                                        o.sat_x, o.sat_y, o.sat_z);
            double tropo = saastamoinenTropoDelay(rough_lat, rough_h, el);

			double iono = 0.0;
			if (!o.used_iono_free){ // Skip Klobuchar for already iono free obse
				iono = klobucharIonoDelay(kp, rough_lat, rough_lon, az, el, epoch.t_gps);
				n_single_freq_epoch++;                        // CORRECTED: moved from Step 1
			} else {
				n_iono_free_epoch++;                          // CORRECTED: moved from Step 1
			}

            Observation o_new = o;
            o_new.pseudorange -= (iono + tropo);
            corrected.push_back(o_new);
        }


        // ── STEP 5: Final solve ──
        PositionSolution sol = solveGpsPosition(corrected, false);
        if (!sol.converged) { failed++; continue; }

		// NEW: compute 3D error vs truth for this epoch
		double err_x = sol.x - TRUE_X;
		double err_y = sol.y - TRUE_Y;
		double err_z = sol.z - TRUE_Z;
		double err_3d = sqrt(err_x*err_x + err_y*err_y + err_z*err_z);
		errors_3d.push_back(err_3d);
        // ── STEP 6: Save to CSV ──
        double lat, lon, h;
        ecefToLla(sol.x, sol.y, sol.z, lat, lon, h);
        // CORRECTED: lat/lon with 9 decimals (was 6). 1e-6 deg = 0.111 m north and
        //   0.111*cos(33.58 deg) = 0.093 m east, a grid visible once errors are ~1 m.
        //   1e-9 deg = 0.1 mm. Every other column keeps 6 decimals, as before.
        csv << fixed << setprecision(6)
            << epoch.t_gps << ","
            << setprecision(9) << lat << "," << lon << ","   // CORRECTED: 9 decimals
            << setprecision(6)                                  // CORRECTED: back to 6
            << h << "," << sol.sigma0 << "," << sol.dop.PDOP << ","
            << corrected.size() << "," << sol.converged << ","
			<< n_iono_free_epoch << "," << n_single_freq_epoch << ","
			<< err_3d << endl;   // NEW columns

        successful++;

	}

	csv.close();

	cout << endl << "=================================" << endl;
	cout << "       PROCESSING COMPLETE" << endl;
	cout << "==========================================" << endl;
	cout << "Total epochs:      " << epochs.size() << endl;
	cout << "Successful fixes:  " << successful << endl;
	cout << "Failed:            " << failed << endl;
	cout << "Success rate:      "
		 << (100.0 * successful / epochs.size()) << "%" << endl;
	cout << "Output: trajectory_dualfreq_crctd.csv" << endl;   // CORRECTED: new file name

	// NEW: 3D error statistics
	if (!errors_3d.empty()) {
		double sum = 0.0, sum_sq = 0.0;
		double min_err = errors_3d[0], max_err = errors_3d[0];
		for (double e : errors_3d) {
			sum += e;
			sum_sq += e * e;
			if (e < min_err) min_err = e;
			if (e > max_err) max_err = e;
		}
		double mean_err = sum / errors_3d.size();
		double rms_err  = sqrt(sum_sq / errors_3d.size());

		cout << endl << "=== 3D ACCURACY vs TRUE POSITION ===" << endl;
		cout << fixed << setprecision(3);
		cout << "Mean 3D error:  " << mean_err << " m" << endl;
		cout << "RMS  3D error:  " << rms_err  << " m" << endl;
		cout << "Min  3D error:  " << min_err  << " m" << endl;
		cout << "Max  3D error:  " << max_err  << " m" << endl;
	}

    return 0;

}
