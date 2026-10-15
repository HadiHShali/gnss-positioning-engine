// real_data_position_withDOP_crctd.cpp - Compute receiver position from real RINEX data,
//                                         with DOP output and an elevation mask
//
// CORRECTED copy of real_data_position_withDOP.cpp (Month 2 / Week 7 / Day 3 review, 2026-10-15).
// The original file stays untouched, so `diff` between the two shows exactly what
// changed. Every change is marked with a "CORRECTED" comment.
//
// THE PROBLEM (time of calculation) - the same bug as Week 6 Day 4:
//   The receiver tags the epoch with the RECEPTION time t_obs (epoch.t_gps). The signal
//   left the satellite about 67-86 ms earlier, at the TRANSMIT time t_tx. In the
//   Earth-fixed frame a GPS satellite moves ~3 km/s, so in that time it travels ~200-300 m.
//   The original evaluated everything at t_obs, so each satellite was placed where it
//   was when the signal ARRIVED, not where the signal LEFT from.
//   Three lines depended on that time:
//     1. Satellite position : computeSatPosECEF(eph, t_obs)            -> now at t_tx
//     2. Satellite clock    : af0 + af1*(t_obs - toe)                   -> now at t_tx,
//                             referenced to toc, + af2 + relativistic term, week wrap
//     3. Earth rotation     : travel time = P/c (contains BOTH clocks)  -> now the true
//                             travel time  (t_obs - dt_rx) - t_tx
//
//   Size of the error (measured in Week 6 Day 4 with the same math):
//     Synthetic BILL data, exact light-time truth: 24.5 m 3D error -> <= 1 mm.
//       - satellite position at the wrong time alone : 22.9 m
//       - missing relativistic clock term alone      :  9.3 m
//       - P/c used as the Sagnac travel time         :  0.4 m per ms of receiver clock
//     Real BILL data: horizontal error 34.4 m -> 6.9 m, RMS residual 28.4 m -> 7.3 m.
//
// HOW THE FIX FITS DAY 3'S TWO PASSES:
//   Original : build obs -> PASS 1 (all sats) -> elevation mask -> PASS 2 (kept sats)
//   Corrected: t_tx, dt_sat and the satellite position are computed ONCE per satellite
//              (they need no receiver information, see computeSatGeom)
//              -> build obs with dt_rx = 0            (receiver clock still unknown)
//              -> PASS 1 (all sats): rough position AND the receiver clock dt_rx
//              -> REBUILD obs with PASS 1's dt_rx     (Sagnac now uses the true flight time)
//              -> elevation mask -> PASS 2 (kept sats): the final, fully corrected solution
//   PASS 1 is only slightly off (Sagnac ~0.4 m per ms of receiver clock). That is far
//   more than good enough for the 10-degree mask and for estimating dt_rx.
//
// NOT CHANGED: satellites.cpp, rinex_obs.cpp, gps_ls_solver.cpp, utils.cpp.
//   computeSatPosECEF() returns the satellite at whatever time it is given, so only the
//   CALLER (this file) passed the wrong time. The helper functions below are copied
//   unchanged from Week 6 Day 4 real_data_position_crctd.cpp (the same math as the
//   tested Week 16 Day 3 sat_geometry.cpp).
//
// ADDED (not a timing issue, marked "ADDED (TGD)"): the L1 group delay TGD.
//   The broadcast clock (af0, af1, af2) is defined for the L1/L2 ionosphere-free
//   combination. A single-frequency L1 code user (this app uses C1C) must shift it to L1:
//       dt_sat(L1) = dt_sat - TGD          (IS-GPS-200, 20.3.3.3.3.2)
//   TGD is a few ns per satellite, i.e. up to ~1-3 m of range. It was parsed into
//   eph.TGD by satellites.cpp but never used. (Dual-frequency ionosphere-free code,
//   Month 3, does NOT need it: there it cancels by design.)
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

// ── CORRECTED: transmit-time helpers (copied unchanged from Week 6 Day 4 _crctd) ──
// Everything inside "namespace { }" is private to this file, like 'static'.
namespace {

    const double MU_GPS = 3.986005e14;       // Earth GM used by GPS (IS-GPS-200) [m^3/s^2]; same value as satellites.cpp
    const double F_REL = -4.442807633e-10;  // relativistic clock constant -2*sqrt(MU)/c^2 [s/sqrt(m)] (IS-GPS-200)
    const double HALF_WEEK = 302400.0;          // half a GPS week [s]
    const double WEEK = 604800.0;          // one GPS week [s]

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
    //   Case                      t                    toc (time of clock)  naive t-toc   wrapped    meaning
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
    //   which wraps tk = t - toe the same way.
    //
    // WHEN IT TRIGGERS: only for epochs near Saturday/Sunday midnight. For BILL at
    //   Sunday noon (t = 43200 s, toc a few hours away) it never changes anything.
    double wrapWeek(double dt)
    {
        if (dt > HALF_WEEK) dt -= WEEK;   // t at week end, toc in next week:  604000 - 0      -> -800 s
        if (dt < -HALF_WEEK) dt += WEEK;   // t at week start, toc in last week:     0 - 597600 -> +7200 s
        return dt;
    }

    // Eccentric anomaly E_k at time t [rad]. Needed only for the relativistic clock term.
    // Same first steps as computeSatPosECEF(): mean motion, mean anomaly, Kepler's equation.
    double eccentricAnomaly(const GpsEphemeris& eph, double t)
    {
        double a = eph.sqrt_a * eph.sqrt_a;                     // semi-major axis [m]
        double n = sqrt(MU_GPS / (a * a * a)) + eph.delta_n;    // corrected mean motion [rad/s]
        double tk = wrapWeek(t - eph.toe);                       // time since toe [s]
        double M = eph.M0 + n * tk;                             // mean anomaly [rad]
        double E = M;                                           // start value (e is small, so E ~ M)
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
        long dow = ((days % 7) + 7) % 7;                        // keep it in 0..6

        return dow * 86400.0 + eph.hour * 3600.0 + eph.minute * 60.0 + eph.second;
    }

    // Everything this app needs about ONE satellite at ONE epoch.
    struct SatGeom {
        double t_tx = 0.0;   // signal transmit time, GPS seconds of week [s]
        double dt_sat = 0.0;   // satellite clock offset at t_tx [s] = polynomial + relativistic (TGD not included)
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
            double t = t_raw - dt_sat;                        // current transmit-time estimate [s]
            double dtc = wrapWeek(t - toc);                     // time since toc [s], week-wrapped
            double poly = eph.clk_bias                           // af0 [s]
                + eph.clk_drift * dtc               // af1 [s/s]
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

    // STEP 2: Parse RINEX OBS (pseudoranges) Noon GPS time
    // (RINEX epochs are in GPS time, so 43200 s of week = Sunday 12:00 GPS, not UTC.)
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
    // receiver's position or clock (the receiver clock cancels, see computeSatGeom),
    // so they stay the same for PASS 1 and PASS 2.
    struct SatInput {
        int prn;        // satellite PRN
        double P;       // measured pseudorange [m]
        SatGeom g;      // transmit time, satellite clock, satellite position at t_tx
        double tgd;     // ADDED (TGD): L1/L2 group delay from the nav record [s]
    };
    vector<SatInput> sats;
    cout << endl << "Building observations..." << endl;

    for (const auto& rec : epoch.gps_records) {
        // Picking the ephemeris with t_obs is fine: 75 ms never changes which toe is nearest.
        const GpsEphemeris* eph = findBestEphemeris(ephs, rec.prn, epoch.t_gps);
        if (!eph) continue;

        // CORRECTED: was  SatPosition sp = computeSatPosECEF(*eph, epoch.t_gps);  (RECEPTION time)
        // CORRECTED: was  sat_clock_correction = af0 + af1*(epoch.t_gps - toe);   (RECEPTION time,
        //            wrong reference epoch, no af2, no relativistic term, no week wrap)
        // Both now come from computeSatGeom(), evaluated at the TRANSMIT time.
        SatInput s;
        s.prn = rec.prn;
        s.P = rec.pseudorange;
        s.g = computeSatGeom(*eph, epoch.t_gps, rec.pseudorange);
        s.tgd = eph->TGD;   // ADDED (TGD): keep this satellite's group delay for buildObservations
        sats.push_back(s);
    }

    // CORRECTED: building the Observation vector is now a small function of the receiver
    // clock dt_rx, because the Earth-rotation (Sagnac) correction needs the TRUE travel time:
    //   travel time = t_rx_true - t_tx = (t_obs - dt_rx) - t_tx
    // The original used P/c, which also contains c*dt_rx and -c*dt_sat. Each millisecond
    // of clock error rotates the satellite ~2 m too far (OMEGA_EARTH * 1e-3 s * 26,560 km),
    // which moved the position ~0.4 m per ms in the Week 6 test.
    // dt_rx is unknown until we solve, so this is called twice: dt_rx = 0 for PASS 1,
    // then PASS 1's clock estimate for PASS 2.
    // [&] means "use the variables of main (sats, epoch) by reference inside this function".
    auto buildObservations = [&](double dt_rx) {
        vector<Observation> out;
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
            // P = rho + c*dt_rx - c*dt_sat(L1)  ->  P + c*dt_sat(L1) = rho + c*dt_rx  (the solver's model)
            // CORRECTED: dt_sat now at t_tx.
            // ADDED (TGD): L1 single-frequency user, so dt_sat(L1) = dt_sat - TGD (IS-GPS-200 20.3.3.3.3.2).
            o.pseudorange = s.P + C_LIGHT * (s.g.dt_sat - s.tgd);
            out.push_back(o);
        }
        return out;
        };

    // CORRECTED: PASS 1 does not know the receiver clock yet, so it uses dt_rx = 0.
    vector<Observation> obs = buildObservations(0.0);
    cout << "Built " << obs.size() << " observation tuples" << endl << endl;

    if (obs.size() < 4) { cerr << "Not enough!" << endl; return 1; }

    // STEP 4: Solve!
    // Pass 1: Solve with all satellites
    // (Rough solution: Sagnac still assumes dt_rx = 0. Used only for the elevation mask
    //  and to estimate dt_rx, both of which it does to far better than needed.)
    PositionSolution sol1 = solveGpsPosition(obs, false);
    cout << "PASS 1: " << obs.size() << " sats, PDOP = " << sol1.dop.PDOP << ", sigma0 = " << sol1.sigma0 << " m" << endl;

    // CORRECTED: rebuild the observations with PASS 1's receiver clock, so the Sagnac
    // rotation now uses the true flight time. Satellite times and clocks don't change;
    // only the rotation angle does (by OMEGA_EARTH * dt_rx).
    obs = buildObservations(sol1.clock_bias);

    // Filter by elevation using rough position
    // (Elevation is computed from the corrected satellite positions and the PASS 1 position.)
    vector<Observation> filtered_obs;
    const double ELEV_MASK = 10.0;  // degrees

    for (const auto& o : obs) {
        double elev = computeElevation(sol1.x, sol1.y, sol1.z,
            o.sat_x, o.sat_y, o.sat_z);
        cout << "  G" << o.prn << "  elevation = " << elev << " deg";
        if (elev >= ELEV_MASK) {
            filtered_obs.push_back(o);
            cout << " KEPT" << endl;
        }
        else {
            cout << " REJECTED (below " << ELEV_MASK << " deg)" << endl;
        }
    }

    // PASS 2: Solve with filtered satellites
    // This is the final solution: transmit-time satellites, satellite clock at t_tx,
    // and Sagnac with the true flight time.
    if (filtered_obs.size() < 4) {
        cerr << "Not Enough Satellites after masking!!" << endl;
        return 1;
    }

    PositionSolution sol2 = solveGpsPosition(filtered_obs, false);
    cout << endl << "PASS 2: " << filtered_obs.size()
        << " sats, PDOP = " << sol2.dop.PDOP
        << ", sigma0 = " << sol2.sigma0 << " m" << endl;

    return 0;
}