// real_data_position_crctd.cpp - Compute receiver position from real RINEX data
//
// CORRECTED copy of real_data_position.cpp (Month 2 / Week 6 review, 2026-10-12).
// The original file stays untouched, so `diff` between the two shows exactly what
// changed. Every change is marked with a "CORRECTED" comment.
//
// THE PROBLEM (time of calculation):
//   The receiver tags the epoch with the RECEPTION time t_obs. The signal left the
//   satellite about 67-86 ms earlier, at the TRANSMIT time t_tx. In the Earth-fixed
//   frame a GPS satellite moves ~3 km/s, so in that time it travels ~200-300 m.
//   The original evaluated everything at t_obs, so each satellite was placed where
//   it was when the signal ARRIVED, not where the signal LEFT from.
//   Three lines depended on that time:
//     1. Satellite position : computeSatPosECEF(eph, t_obs)          -> now at t_tx
//     2. Satellite clock    : af0 + af1*(t_obs - toe)                 -> now at t_tx,
//                             referenced to toc, + af2 + relativistic term, week wrap
//     3. Earth rotation     : travel time = P/c (contains BOTH clocks)  -> now the true
//                             travel time  (t_obs - dt_rx) - t_tx
//
// NOT CHANGED: satellites.cpp, rinex_obs.cpp, gps_ls_solver.cpp, utils.cpp.
//   They are correct for this problem: computeSatPosECEF() returns the satellite at
//   whatever time it is given, so only the CALLER (this file) passed the wrong time.
//   The math below is the same as the tested Week 16 Day 3 helper (sat_geometry.cpp).
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "gps_ls_solver.h"        // Week 6 solver
#include "utils.h"                  // ecefToLla
#include "satellites.h"            // Week 5 — parseRinexNav, computeSatPosECEF
#include "rinex_obs.h"             // NEW today
const double C_LIGHT = 299792458.0;
const double OMEGA_EARTH = 7.2921151467e-5;
using namespace std;

// ── CORRECTED: transmit-time helpers (same math as Week 16 sat_geometry.cpp) ──
// Everything inside "namespace { }" is private to this file, like 'static'.
namespace {

const double MU_GPS    = 3.986005e14;       // Earth GM used by GPS (IS-GPS-200) [m^3/s^2]; same value as satellites.cpp
const double F_REL     = -4.442807633e-10;  // relativistic clock constant -2*sqrt(MU)/c^2 [s/sqrt(m)] (IS-GPS-200)
const double HALF_WEEK = 302400.0;          // half a GPS week [s]
const double WEEK      = 604800.0;          // one GPS week [s]

// Wrap a time DIFFERENCE into [-half week, +half week].
// Times are "seconds of week", so at the start of a week t is ~0 while toe/toc
// may still belong to the previous week (~597600 s). Without the wrap,
// t - toc = -597600 s and the clock drift term would be multiplied by a week.
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
// toc is the date/time on the FIRST line of the nav record (eph.year ... eph.second).
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
    double dt_sat = 0.0;   // satellite clock offset at t_tx [s] = polynomial + relativistic (TGD not included)
    double dt_rel = 0.0;   // relativistic part of dt_sat [s] (kept separately only to print it)
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

} // namespace
// ─────────────────────────────────────────────────────────────────────────────

int main() {
    cout << "=== REAL GPS POSITIONING (transmit-time corrected) ===" << endl << endl;   // CORRECTED: title only

    // STEP 1: Parse RINEX NAV (broadcast ephemerides)
    string nav_file = "../../data/BRDC00IGS_R_20240070000_01D_MN.rnx";
    vector<GpsEphemeris> ephs = parseRinexNav(nav_file);
    cout << "Parsed " << ephs.size() << " ephemerides" << endl;

    // STEP 2: Parse RINEX OBS (pseudoranges) Noon utc
    string obs_file = "../../data/BILL00USA_R_20240070000_01D_15S_MO.rnx";
    ObsEpoch epoch = parseRinexObsOneEpoch(obs_file, 43200.0);
    cout << "Epoch: " << epoch.year << "-" << epoch.month << "-" << epoch.day
        << " " << epoch.hour << ":" << epoch.minute << ":" << epoch.second << endl;
    cout << "GPS sec of week: " << epoch.t_gps << endl;
    cout << "Number of GPS observations: " << epoch.gps_records.size() << endl;

    if (epoch.gps_records.size() < 4) {
        cerr << "Need at least 4 satellites!" << endl;
        return 1;
    }


    // STEP 3: CORRECTED — for each observed satellite, compute its state at the TRANSMIT time.
    // Done ONCE per satellite: t_tx, dt_sat and the position do not depend on the
    // receiver's position or clock (the receiver clock cancels, see computeSatGeom).
    struct SatInput {
        int prn;                        // satellite PRN
        double P;                       // measured pseudorange [m]
        const GpsEphemeris* eph;        // ephemeris used (kept to print the timing table)
        SatGeom g;                      // transmit time, clock, position at t_tx
    };
    vector<SatInput> sats;
    cout << endl << "Building observations..." << endl;

    for (const auto& rec : epoch.gps_records) {
        // Picking the ephemeris with t_obs is fine: 75 ms never changes which toe is nearest.
        const GpsEphemeris* eph = findBestEphemeris(ephs, rec.prn, epoch.t_gps);
        if (!eph) continue;

        // CORRECTED: was computeSatPosECEF(*eph, epoch.t_gps), i.e. at the RECEPTION time.
        // CORRECTED: the satellite clock is now also evaluated at t_tx, from toc, with af2,
        //            the relativistic term and the week wrap (was af0 + af1*(t_obs - toe)).
        SatGeom g = computeSatGeom(*eph, epoch.t_gps, rec.pseudorange);
        sats.push_back({ rec.prn, rec.pseudorange, eph, g });
    }
    if (sats.size() < 4) { cerr << "Not enough!" << endl; return 1; }

    // STEP 4: CORRECTED — Earth rotation (Sagnac) needs the TRUE travel time, and solve.
    //   travel time = t_rx_true - t_tx = (t_obs - dt_rx) - t_tx
    // The original used P/c, which also contains c*dt_rx and -c*dt_sat. Each millisecond
    // of clock error rotates the satellite ~2 m too far (OMEGA_EARTH * 1e-3 * 26,560 km).
    // dt_rx is unknown until we solve, so solve twice:
    //   pass 1: assume dt_rx = 0  -> Sagnac off by ~0.4 m of position per ms of dt_rx; dt_rx estimated
    //   pass 2: use the pass-1 dt_rx -> travel time now exact to ~ns: final solution
    double dt_rx = 0.0;                                      // receiver clock bias [s], unknown at first
    PositionSolution sol = {};
    double x1 = 0.0, y1 = 0.0, z1 = 0.0;                     // pass-1 position, to show the pass-2 change
    vector<Observation> obs;

    for (int pass = 1; pass <= 2; ++pass) {
        obs.clear();
        for (const auto& s : sats) {
            // ── EARTH ROTATION CORRECTION (Sagnac) ──
            // CORRECTED: was  travel_time = rec.pseudorange / C_LIGHT
            double travel_time = (epoch.t_gps - dt_rx) - s.g.t_tx;   // true signal flight time [s]
            double rotation_angle = OMEGA_EARTH * travel_time;       // Earth's rotation during the flight [rad]
            double cos_r = cos(rotation_angle);
            double sin_r = sin(rotation_angle);
            // Rotate the satellite from the Earth-fixed frame of t_tx into the frame of
            // the reception time (same rotation as the original, only the angle changed).
            double sx_new = cos_r * s.g.X + sin_r * s.g.Y;
            double sy_new = -sin_r * s.g.X + cos_r * s.g.Y;
            // ───────────────────────────────────────

            Observation o;
            o.prn = s.prn;
            o.sat_x = sx_new;
            o.sat_y = sy_new;
            o.sat_z = s.g.Z;                                         // rotation about Z leaves Z unchanged
            // P = rho + c*dt_rx - c*dt_sat  ->  P + c*dt_sat = rho + c*dt_rx  (the solver's model)
            o.pseudorange = s.P + C_LIGHT * s.g.dt_sat;              // CORRECTED: dt_sat now at t_tx
            obs.push_back(o);
        }
        if (pass == 1) cout << "Built " << obs.size() << " observation tuples" << endl << endl;

        sol = solveGpsPosition(obs, false);
        if (!sol.converged) {
            cerr << "Did not converge! (pass " << pass << ")" << endl;
            return 1;
        }
        dt_rx = sol.clock_bias;                                      // receiver clock for the next pass
        if (pass == 1) { x1 = sol.x; y1 = sol.y; z1 = sol.z; }
    }

    // CORRECTED: per-satellite timing table, so the size of each correction is visible.
    //   flight  : true signal travel time
    //   moved   : how far the satellite is from where the ORIGINAL put it (at t_obs)
    //   dt_sat  : satellite clock at t_tx;  rel: its relativistic part, in meters
    cout << "PRN | flight [ms] | moved [m] | dt_sat [us] | rel [m]" << endl;
    cout << "----+-------------+-----------+-------------+--------" << endl;
    for (const auto& s : sats) {
        SatPosition p_rx = computeSatPosECEF(*s.eph, epoch.t_gps);  // where the ORIGINAL placed it
        double mx = p_rx.X - s.g.X, my = p_rx.Y - s.g.Y, mz = p_rx.Z - s.g.Z;
        cout << "G" << setfill('0') << setw(2) << s.prn << setfill(' ')
            << " | " << fixed << setprecision(3) << setw(11) << ((epoch.t_gps - dt_rx) - s.g.t_tx) * 1e3
            << " | " << setprecision(1) << setw(9) << sqrt(mx * mx + my * my + mz * mz)
            << " | " << setprecision(3) << setw(11) << s.g.dt_sat * 1e6
            << " | " << setw(6) << C_LIGHT * s.g.dt_rel << endl;
    }
    double d1 = sqrt((sol.x - x1) * (sol.x - x1) + (sol.y - y1) * (sol.y - y1) + (sol.z - z1) * (sol.z - z1));
    cout << endl << "Pass 1 -> pass 2 position change (receiver clock in travel time): "
        << setprecision(3) << d1 << " m" << endl << endl;

    // STEP 5: Print results
    cout << "=== POSITION SOLUTION ===" << endl;
    cout << fixed << setprecision(3);
    cout << "Iterations: " << sol.iterations << endl;
    cout << "ECEF:" << endl;
    cout << "  X = " << sol.x << " m" << endl;
    cout << "  Y = " << sol.y << " m" << endl;
    cout << "  Z = " << sol.z << " m" << endl;

    double lat, lon, h;
    ecefToLla(sol.x, sol.y, sol.z, lat, lon, h);
    cout << endl << "LLA:" << endl;
    cout << setprecision(6);
    cout << "  Latitude:  " << lat << " deg" << endl;
    cout << "  Longitude: " << lon << " deg" << endl;
    cout << setprecision(2);
    cout << "  Height:    " << h << " m" << endl;

    cout << endl << "Clock bias: " << sol.clock_bias * 1e6 << " microseconds" << endl;
    cout << "RMS residual: " << sol.rms_residual << " m" << endl;
    cout << "Sigma-0:      " << sol.sigma0 << " m" << endl;

    return 0;
}