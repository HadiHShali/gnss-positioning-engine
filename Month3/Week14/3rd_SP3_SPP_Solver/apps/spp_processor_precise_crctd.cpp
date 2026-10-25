// spp_processor_precise_crctd.cpp - Dual-frequency SPP with IGS precise orbits (SP3) and clocks (CLK)
//
// CORRECTED copy of spp_processor_precise.cpp (Month 3 / Week 14 review, 2026-10-24).
// The original file stays untouched, so `diff` between the two shows exactly what
// changed. Every change is marked with "CORRECTED".
//
// WHAT THE ORIGINAL ALREADY HAD RIGHT:
//   - precise orbit + precise clock per satellite, broadcast as fallback
//   - relativistic periodic term added to the PRECISE clock too (IGS clocks do not
//     contain it - your Week 14 finding)
//   - no TGD for P_IF: IGS clocks, like the broadcast clock, refer to the code
//     ionosphere-free combination
//
// WHAT WAS STILL WRONG:
//   1-3. The timing bug of Weeks 6-13, now in the SP3/CLK calls:
//        orbit interpolated at the RECEPTION time t_obs   -> now at the transmit time t_tx
//        clock interpolated at t_obs                       -> now at t_tx
//        Sagnac travel time = P_IF/c (contains both clocks) -> true travel time (t_obs - dt_rx) - t_tx
//        The broadcast fallback gets the same fix (computeSatGeom, unchanged from Week 13).
//   4. TGD was skipped for SINGLE-frequency satellites when the precise clock was used
//      ("NO TGD - precise clock referenced to iono-free"). That reason is exactly why
//      TGD IS needed for P1 alone: the IGS clock is the IF clock, and TGD converts it to
//      the L1 P(Y) clock, same as with the broadcast clock. Affects only satellites
//      without C2W (n_single_freq); c*TGD is up to a few metres.
//   5. Mixed sources: a precise POSITION could be paired with the BROADCAST clock when
//      the CLK file had no value. Orbit and clock errors of one product cancel partly in
//      the range; mixing products breaks that (~1 m). Now a satellite uses precise
//      orbit AND clock, or broadcast orbit AND clock.
//   6. Week 13 output fixes: counters after the mask, 9-decimal lat/lon, new file name.
//   7. Includes interpolation_crctd.h: the original interpolation.h extrapolated the orbit
//      past the last SP3 epoch (23:45, up to 1.6 m error by 24:00) and interpolated across
//      data gaps. Those satellite-epochs now fall back to broadcast orbit AND clock.
//      (Using the NEXT day's SP3/CLK files would keep them precise - a PPP-week item.)
//
// KNOWN, NOT CHANGED (needs an ANTEX file, planned for the PPP weeks):
//   SP3 orbits and CLK clocks refer to the satellite CENTRE OF MASS; broadcast orbits
//   refer to the antenna PHASE CENTRE. The satellite antenna offset (~1-2 m, mostly
//   along the satellite-to-Earth direction, different per satellite block) is not
//   applied. Most of it is common to all satellites and goes into the receiver
//   clock; the rest is cm to a few dm.
//
// NOT CHANGED: satellites.cpp, rinex_obs.cpp, gps_ls_solver.cpp, utils.cpp,
//   sp3_reader.cpp, clk_reader.cpp (all latent-only issues, see the review notes).
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include "gps_ls_solver.h"
#include "utils.h"
#include "satellites.h"
#include "rinex_obs.h"
#include "sp3_reader.h"
#include "clk_reader.h"
#include "interpolation_crctd.h"   // CORRECTED: was "interpolation.h" (no extrapolation past 23:45, gap checks)

using namespace std;
using namespace Eigen;

const double C_LIGHT = 299792458.0;
const double OMEGA_EARTH = 7.2921151467e-5;

// ── CORRECTED: transmit-time helpers (copied unchanged from Week 13 _crctd = Week 9) ──
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
//   which wraps tk = t - toe the same way (and in the original's clock block, lines 117-118).
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

// CORRECTED: the precise-product version of computeSatGeom(). Same three steps,
// same transmit-time logic; only the SOURCE of the clock and the position changes:
//   clock    : CLK file (IGS, refers to the P1/P2 ionosphere-free combination)
//              + relativistic term (IGS clocks do not contain it), at the transmit time
//   position : SP3 file (Lagrange interpolation), at the transmit time
// The relativistic term still uses the broadcast orbit (eccentricity and E): it only
// needs e and E, and the broadcast values are far more accurate than it requires.
// Returns false if the CLK or SP3 file has no value for this satellite at t_tx; the
// caller then uses the broadcast orbit AND clock, never one of each.
bool computeSatGeomPrecise(const vector<Sp3Epoch>& sp3, const vector<ClkEpoch>& clk,
                           const string& prn, const GpsEphemeris& eph,
                           double t_obs, double P, SatGeom& g)
{
    // 1. Pseudorange-based transmit time (receiver clock cancels, see computeSatGeom)
    const double t_raw = t_obs - P / C_LIGHT;

    // 2. Satellite clock at the transmit time; iterate twice because t_tx needs dt_sat
    double dt_sat = 0.0, dt_rel = 0.0;
    for (int k = 0; k < 2; ++k) {
        double t = t_raw - dt_sat;                           // current transmit-time estimate [s]
        double clk_t;                                        // IGS clock at t [s]
        if (!interpolateSatClock(clk, prn, t, clk_t)) return false;
        dt_rel = F_REL * eph.e * eph.sqrt_a                  // relativistic term, same as broadcast
               * sin(eccentricAnomaly(eph, t));
        dt_sat = clk_t + dt_rel;
    }
    g.dt_sat = dt_sat;
    g.dt_rel = dt_rel;

    // 3. True transmit time, and the precise position AT THAT TIME
    g.t_tx = t_raw - dt_sat;
    double pos3[3];
    if (!interpolateSatPosition(sp3, prn, g.t_tx, pos3)) return false;
    g.X = pos3[0];  g.Y = pos3[1];  g.Z = pos3[2];           // centre of mass (see header note)
    return true;
}

// CORRECTED: one satellite's inputs for one epoch, computed once and reused by both solves.
// Same as Week 13 _crctd.
struct SatInput {
    int prn;          // satellite PRN
    double P;         // pseudorange used in the solve [m]: P_IF (dual) or P1 = C1W (single)
    bool iono_free;   // true -> P is P_IF: no TGD, no Klobuchar
    SatGeom g;        // transmit time, satellite clock, satellite position at t_tx
    double tgd;       // L1/L2 group delay from the nav record [s]
};

// CORRECTED: build the solver's Observation vector for one epoch, for a given receiver clock.
// Identical to Week 13 _crctd: true-travel-time Sagnac, TGD only for single frequency.
// The TGD rule is the same for precise and broadcast clocks, since BOTH refer to the
// P1/P2 ionosphere-free combination (CORRECTED item 4).
vector<Observation> buildObservations(const vector<SatInput>& sats, double t_obs, double dt_rx)
{
    vector<Observation> out;
    for (const auto& s : sats) {
        // Sagnac
        // CORRECTED: was  travel_time = base_pr / C_LIGHT
        double travel_time = (t_obs - dt_rx) - s.g.t_tx;            // true signal flight time [s]
        double rot = OMEGA_EARTH * travel_time;                     // Earth's rotation during the flight [rad]
        double cr = cos(rot), sr = sin(rot);
        double sx =  cr * s.g.X + sr * s.g.Y;
        double sy = -sr * s.g.X + cr * s.g.Y;

        Observation o;
        o.prn = s.prn;
        o.sat_x = sx; o.sat_y = sy; o.sat_z = s.g.Z;                // rotation about Z leaves Z unchanged
        //   P_IF : dt_sat(IF) = dt_sat           (precise or broadcast clock IS the IF clock)
        //   P1   : dt_sat(L1) = dt_sat - TGD     (CORRECTED for precise clocks, item 4)
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
	cout << "=== FULL SPP PROCESSOR (precise SP3/CLK, transmit-time corrected) ===" << endl;   // CORRECTED: title only

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


		    // NEW: parse precise products
    auto sp3_epochs = parseSp3("../../data/orbit.sp3");
    auto clk_epochs = parseClk("../../data/clock.clk");
    cout << "Parsed " << sp3_epochs.size() << " SP3 epochs, "
         << clk_epochs.size() << " CLK epochs" << endl;


	// Open trajectoory CSV
	// CORRECTED (output file name only): was "../../data/trajectory_dualfreq_SP3CLK.csv".
	ofstream csv("../../data/trajectory_dualfreq_SP3CLK_crctd.csv");
	csv << "t_gps,lat_deg,lon_deg,height_m,sigma0_m,pdop,n_sats,converged,n_iono_free,n_single_freq,err_3d" << endl;

	// Track statistics
	int successful = 0;
	int failed = 0;
	int n_precise_pos = 0, n_broadcast_pos = 0;   // satellite-epochs by orbit/clock source
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

		// Step1: CORRECTED — satellite state at the TRANSMIT time, once per satellite,
		// from the precise products when both orbit and clock exist, else from broadcast.
		vector<SatInput> sats;
		for (const auto& rec: epoch.gps_records)
		{
			// Still needed: TGD, the relativistic term (e, E) and the fallback.
			// Picking it with t_obs is fine: 75 ms never changes which toe is nearest.
			const GpsEphemeris* eph = findBestEphemeris(ephs, rec.prn, epoch.t_gps);
			if (!eph) continue;

			string prn = "G" + string(rec.prn < 10 ? "0" : "") + to_string(rec.prn);

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
			// CORRECTED: the per-epoch iono-free/single counters moved to Step 4 (after the mask).

			// CORRECTED: was  interpolateSatPosition(sp3_epochs, prn, epoch.t_gps, ...)   (RECEPTION time)
			//            and  interpolateSatClock(clk_epochs, prn, epoch.t_gps, ...)       (RECEPTION time)
			//            and  the broadcast fallback at epoch.t_gps from toe, and a precise
			//            orbit could be paired with the broadcast clock.
			// Now: precise orbit AND clock at t_tx, or broadcast orbit AND clock at t_tx.
			// The transmit time comes from the measured P1 (C1W), as in Week 13.
			SatInput s;
			s.prn       = rec.prn;
			s.P         = base_pr;                 // P_IF or P1: what the solve uses
			s.iono_free = is_iono_free;
			s.tgd       = eph->TGD;                // applied only when !iono_free (buildObservations)
			if (computeSatGeomPrecise(sp3_epochs, clk_epochs, prn, *eph,
			                          epoch.t_gps, rec.pseudorange, s.g)) {
				n_precise_pos++;
			} else {
				s.g = computeSatGeom(*eph, epoch.t_gps, rec.pseudorange);   // broadcast fallback
				n_broadcast_pos++;
			}
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

		// CORRECTED: rebuild with the rough solve's receiver clock (true-travel-time Sagnac).
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
        // CORRECTED: lat/lon with 9 decimals (was 6, a 0.11 m x 0.09 m grid). Same as Week 13.
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
		 // in the final summary block:
	cout << "Precise orbit+clock used: " << n_precise_pos << endl;     // CORRECTED: label (both, not orbit only)
	cout << "Broadcast fallback used:  " << n_broadcast_pos << endl;
	cout << "Output: trajectory_dualfreq_SP3CLK_crctd.csv" << endl;   // CORRECTED: new file name

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
