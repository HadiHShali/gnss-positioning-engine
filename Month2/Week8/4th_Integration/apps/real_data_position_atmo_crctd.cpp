// real_data_position_atmo_crctd.cpp - Compute receiver position from real RINEX data,
//                                     with Klobuchar + Saastamoinen atmospheric corrections
//
// CORRECTED copy of real_data_position_atmo.cpp (Month 2 / Week 8 / last day review, 2026-10-15).
// The original file stays untouched, so `diff` between the two shows exactly what
// changed. Every change is marked with "CORRECTED".
//
// WHAT THE ORIGINAL ALREADY HAD RIGHT (added in Week 8):
//   the satellite clock now had af2, the relativistic term, the week wrap and TGD.
//
// WHAT WAS STILL WRONG (time of calculation) - the same bug as Week 6 Day 4 and Week 7:
//   The receiver tags the epoch with the RECEPTION time t_obs (epoch.t_gps). The signal
//   left the satellite about 67-86 ms earlier, at the TRANSMIT time t_tx. In the
//   Earth-fixed frame a GPS satellite moves ~3 km/s, so in that time it travels ~200-300 m.
//     1. Satellite position : computeSatPosECEF(eph, t_obs)            -> now at t_tx
//                             (the big one: 22.9 m of position error in the Week 6 test)
//     2. Earth rotation     : travel time = P/c (contains BOTH clocks)  -> now the true
//                             travel time (t_obs - dt_rx) - t_tx  (0.4 m per ms of receiver clock)
//     3. Satellite clock    : evaluated at t_obs and referenced to toe  -> now at t_tx, from toc
//                             (numerically tiny here: 75 ms x af1 ~ 1e-12 s ~ 0.3 mm, and toc = toe
//                              for GPS. Changed anyway so all three use the same tested helper.)
//
// HOW THE FIX FITS THIS FILE'S PASSES:
//   Original : build obs -> PASS 1a (all sats) -> mask -> PASS 1 (masked, no atmo) -> atmo -> PASS 2
//   Corrected: t_tx, dt_sat and the satellite position are computed ONCE per satellite
//              (they need no receiver information, see computeSatGeom)
//              -> build obs with dt_rx = 0            (receiver clock still unknown)
//              -> PASS 1a (all sats): rough position AND the receiver clock dt_rx
//              -> REBUILD obs with PASS 1a's dt_rx    (Sagnac now uses the true flight time)
//              -> mask -> PASS 1 -> atmosphere -> PASS 2   (unchanged from here on)
//   The atmosphere shifts dt_rx by only tens of ns, which changes the Sagnac rotation by
//   < 0.1 mm, so one rebuild after PASS 1a is enough.
//
// EXPECTED CHECK: PASS 1 (masked, no atmo) uses the same epoch, the same 8 satellites and the
//   same satellite math (incl. TGD) as Week 7 Day 3 real_data_position_withDOP_crctd.cpp,
//   so its sigma0 should come out close to that run's 0.61 m.
//
// CORRECTED (truth, not a timing issue): TRUE_X/Y/Z are now BILL00USA's IGS coordinates
//   (-2420420.501, -4737131.551, 3507827.665), the ground truth confirmed in Week 6.
//   The original's values were 1.8 m away (dX = -1.54 m, dY = -0.92 m, dZ = -0.06 m).
//
// FLAGGED, NOT CHANGED (not timing issues):
//   - "clock bias = ... m" prints sol1.clock_bias, which the solver stores in SECONDS.
//   - Correct as they are: atmosphere is subtracted (P - I - T), and Klobuchar uses the
//     reception time, which is the right time for the ionospheric pierce point.
//
// NOT CHANGED: satellites.cpp, rinex_obs.cpp, gps_ls_solver.cpp, utils.cpp.
//   The helper functions below are copied unchanged from Week 7 Day 4 _crctd
//   (= Week 6 Day 4 = Week 16 Day 3 sat_geometry.cpp math).
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "gps_ls_solver.h"        // Week 6 solver
#include "utils.h"                 // ecefToLla
#include "satellites.h"            // parseRinexNav, computeSatPosECEF, atmospheric models
#include "rinex_obs.h"             // RINEX observation parser

const double C_LIGHT = 299792458.0;
const double OMEGA_EARTH = 7.2921151467e-5;
using namespace std;

// ── CORRECTED: transmit-time helpers (copied unchanged from Week 7 Day 4 _crctd) ──
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
    //   which wraps tk = t - toe the same way (and in the original's clock block, lines 75-76).
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
        long dow = ((days % 7) + 7) % 7;                        // keep it in 0..6

        return dow * 86400.0 + eph.hour * 3600.0 + eph.minute * 60.0 + eph.second;
    }

    // Everything this app needs about ONE satellite at ONE epoch.
    struct SatGeom {
        double t_tx = 0.0;   // signal transmit time, GPS seconds of week [s]
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
    cout << "=== REAL GPS POSITIONING WITH ATMOSPHERIC CORRECTIONS (transmit-time corrected) ===" << endl << endl;   // CORRECTED: title only

    // ─────────────────────────────────────────────────────────────
    // STEP 0: Parse RINEX NAV (broadcast ephemerides)
    // ─────────────────────────────────────────────────────────────
    string nav_file = "../../data/BRDC00IGS_R_20240070000_01D_MN.rnx";
    vector<GpsEphemeris> ephs = parseRinexNav(nav_file);
    cout << "Parsed " << ephs.size() << " ephemerides" << endl;

    // ─────────────────────────────────────────────────────────────
    // STEP 1: Parse Klobuchar iono parameters
    // ─────────────────────────────────────────────────────────────
    KlobucharParams kp = parseKlobucharFromNav(nav_file);
    if (kp.valid) {
        cout << "Loaded Klobuchar Iono parameters" << endl;
    }
    else {
        cout << "WARNING: no Klobuchar params - iono correction disabled" << endl;
    }

    // ─────────────────────────────────────────────────────────────
    // STEP 2: Parse RINEX OBS (pseudoranges) at noon GPS time
    // (RINEX epochs are in GPS time, so 43200 s of week = Sunday 12:00 GPS, not UTC.)
    // ─────────────────────────────────────────────────────────────
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

    // ─────────────────────────────────────────────────────────────
    // STEP 3: CORRECTED — satellite state at the TRANSMIT time, once per satellite
    // ─────────────────────────────────────────────────────────────
    // t_tx, dt_sat and the satellite position do not depend on the receiver's position
    // or clock (the receiver clock cancels, see computeSatGeom), so they are computed
    // once here and reused by every pass below.
    struct SatInput {
        int prn;        // satellite PRN
        double P;       // measured pseudorange [m]
        SatGeom g;      // transmit time, satellite clock, satellite position at t_tx
        double tgd;     // L1/L2 group delay from the nav record [s] (the original applied it too)
    };
    vector<SatInput> sats;
    cout << endl << "Building observations..." << endl;

    for (const auto& rec : epoch.gps_records) {
        // Picking the ephemeris with t_obs is fine: 75 ms never changes which toe is nearest.
        const GpsEphemeris* eph = findBestEphemeris(ephs, rec.prn, epoch.t_gps);
        if (!eph) continue;

        // CORRECTED: was  SatPosition sp = computeSatPosECEF(*eph, epoch.t_gps);   (RECEPTION time)
        // CORRECTED: the original's clock block (af0 + af1*dt + af2*dt^2 + relativistic,
        //            week-wrapped) was evaluated at the RECEPTION time and referenced to toe.
        //            computeSatGeom() evaluates the same terms at the TRANSMIT time, from toc.
        SatInput s;
        s.prn = rec.prn;
        s.P = rec.pseudorange;
        s.g = computeSatGeom(*eph, epoch.t_gps, rec.pseudorange);
        s.tgd = eph->TGD;   // TGD (Timing Group Delay) - for single-frequency L1 users, as in the original
        sats.push_back(s);
    }

    // CORRECTED: building the Observation vector is now a small function of the receiver
    // clock dt_rx, because the Earth-rotation (Sagnac) correction needs the TRUE travel time:
    //   travel time = t_rx_true - t_tx = (t_obs - dt_rx) - t_tx
    // The original used P/c, which also contains c*dt_rx and -c*dt_sat. Each millisecond
    // of clock error rotates the satellite ~2 m too far (OMEGA_EARTH * 1e-3 s * 26,560 km),
    // which moved the position ~0.4 m per ms in the Week 6 test.
    // dt_rx is unknown until we solve, so this is called twice: dt_rx = 0 for PASS 1a,
    // then PASS 1a's clock estimate for everything after it.
    // [&] means "use the variables of main (sats, epoch) by reference inside this function".
    auto buildObservations = [&](double dt_rx) {
        vector<Observation> out;
        for (const auto& s : sats) {
            // Earth rotation correction (Sagnac)
            // CORRECTED: was  travel_time = rec.pseudorange / C_LIGHT
            double travel_time = (epoch.t_gps - dt_rx) - s.g.t_tx;   // true signal flight time [s]
            double rotation_angle = OMEGA_EARTH * travel_time;       // Earth's rotation during the flight [rad]
            double cos_r = cos(rotation_angle);
            double sin_r = sin(rotation_angle);
            // Rotate the satellite from the Earth-fixed frame of t_tx into the frame of
            // the reception time (same rotation as the original, only the angle changed).
            double sx_new = cos_r * s.g.X + sin_r * s.g.Y;
            double sy_new = -sin_r * s.g.X + cos_r * s.g.Y;

            Observation o;
            o.prn = s.prn;
            o.sat_x = sx_new;
            o.sat_y = sy_new;
            o.sat_z = s.g.Z;                                         // rotation about Z leaves Z unchanged
            // P = rho + c*dt_rx - c*dt_sat(L1)  ->  P + c*dt_sat(L1) = rho + c*dt_rx  (the solver's model)
            // CORRECTED: dt_sat now at t_tx, from toc.
            // TGD as in the original: dt_sat(L1) = dt_sat - TGD (IS-GPS-200 20.3.3.3.3.2).
            o.pseudorange = s.P + C_LIGHT * (s.g.dt_sat - s.tgd);
            out.push_back(o);
        }
        return out;
        };

    // CORRECTED: PASS 1a does not know the receiver clock yet, so it uses dt_rx = 0.
    vector<Observation> obs = buildObservations(0.0);
    cout << "Built " << obs.size() << " observation tuples" << endl << endl;

    if (obs.size() < 4) { cerr << "Not enough observations!" << endl; return 1; }

    // ─────────────────────────────────────────────────────────────
    // PASS 1a: Dirty rough solve with all satellites (to get initial position)
    // (and, CORRECTED, the receiver clock dt_rx for the Sagnac rebuild below)
    // ─────────────────────────────────────────────────────────────
    PositionSolution sol1a = solveGpsPosition(obs, false);
    cout << "PASS 1a (all sats): " << obs.size() << " sats, PDOP = "
        << sol1a.dop.PDOP << ", sigma0 = " << sol1a.sigma0 << " m" << endl;

    // CORRECTED: rebuild the observations with PASS 1a's receiver clock, so the Sagnac
    // rotation now uses the true flight time. Satellite times and clocks don't change;
    // only the rotation angle does (by OMEGA_EARTH * dt_rx). Everything below uses these.
    obs = buildObservations(sol1a.clock_bias);

    // ─────────────────────────────────────────────────────────────
    // Filter by elevation using PASS 1a rough position
    // ─────────────────────────────────────────────────────────────
    vector<Observation> filtered_obs;
    const double ELEV_MASK = 10.0;  // degrees

    for (const auto& o : obs) {
        double elev = computeElevation(sol1a.x, sol1a.y, sol1a.z,
            o.sat_x, o.sat_y, o.sat_z);
        if (elev >= ELEV_MASK) {
            filtered_obs.push_back(o);
        }
        else {
            cout << "  G" << o.prn << " rejected (elev = " << elev << " deg)" << endl;
        }
    }

    if (filtered_obs.size() < 4) {
        cerr << "Not enough satellites after masking!" << endl;
        return 1;
    }

    // ─────────────────────────────────────────────────────────────
    // PASS 1: CLEAN rough solve with masked satellites (no atmo yet)
    // (CORRECTED satellites: expect sigma0 close to Week 7 Day 3 _crctd's 0.61 m)
    // ─────────────────────────────────────────────────────────────
    PositionSolution sol1 = solveGpsPosition(filtered_obs, false);
    cout << endl << "PASS 1 (masked, no atmo): " << filtered_obs.size()
        << " sats, PDOP = " << sol1.dop.PDOP
        << ", sigma0 = " << sol1.sigma0 << " m" << endl;

    // Convert clean rough position to lat/lon/height for atmo models
    double rough_lat, rough_lon, rough_h;
    ecefToLla(sol1.x, sol1.y, sol1.z, rough_lat, rough_lon, rough_h);
    cout << "Rough position lat=" << rough_lat << " lon=" << rough_lon
        << " h=" << rough_h << endl << endl;


    // After computing sol1:
    // CORRECTED (truth): BILL00USA IGS coordinates, the ground truth confirmed in Week 6.
    //   The original used (-2420422.0410, -4737132.4707, 3507827.6034), which is 1.8 m away
    //   (E -0.95, N -0.89, U +1.23 m). Scores from the original and from this app are
    //   therefore NOT against the same truth any more.
    // Use the XYZ, not a listed lat/lon/h: the XYZ converts to h = 470.061 m (WGS84),
    //   while the listed h = 470.478 m does not match it (0.42 m apart), and a lat/lon
    //   rounded to 3 decimals (33.578, -117.065) is only good to ~100 m.
    const double TRUE_X = -2420420.501;
    const double TRUE_Y = -4737131.551;
    const double TRUE_Z = 3507827.665;

    double sol1_err_x = sol1.x - TRUE_X;
    double sol1_err_y = sol1.y - TRUE_Y;
    double sol1_err_z = sol1.z - TRUE_Z;

    cout << "PASS 1 sol vs truth:" << endl;
    cout << "  dX = " << sol1_err_x << " m" << endl;
    cout << "  dY = " << sol1_err_y << " m" << endl;
    cout << "  dZ = " << sol1_err_z << " m" << endl;
    cout << "  clock bias = " << sol1.clock_bias << " m" << endl;   // FLAGGED (not changed): this value is in SECONDS
    cout << "  iterations = " << sol1.iterations << endl;
    cout << "  converged  = " << sol1.converged << endl;

    cout << "PASS 1 residuals:" << endl;
    for (int i = 0; i < sol1.residuals.size(); i++) {
        cout << "  " << i << ": " << sol1.residuals(i) << " m" << endl;
    }


    cout << "PASS 1 ECEF: X=" << sol1.x
        << " Y=" << sol1.y
        << " Z=" << sol1.z << endl;
    cout << "Truth  ECEF: X=" << TRUE_X
        << " Y=" << TRUE_Y
        << " Z=" << TRUE_Z << endl;
    // ─────────────────────────────────────────────────────────────
    // PASS 2: Apply atmospheric corrections to SAME masked satellites
    // (Unchanged: the atmosphere is SUBTRACTED, P - I - T, and Klobuchar uses the
    //  reception time t_obs, which is the right time for the ionospheric pierce point.)
    // ─────────────────────────────────────────────────────────────
    vector<Observation> corrected_obs;
    cout << "Applying atmospheric corrections:" << endl;

    for (const auto& o : filtered_obs) {
        double el = computeElevation(sol1.x, sol1.y, sol1.z,
            o.sat_x, o.sat_y, o.sat_z);
        double az = computeAzimuth(sol1.x, sol1.y, sol1.z,
            o.sat_x, o.sat_y, o.sat_z);

        double iono = klobucharIonoDelay(kp, rough_lat, rough_lon,
            az, el, epoch.t_gps);
        double tropo = saastamoinenTropoDelay(rough_lat, rough_h, el);

        Observation o_new = o;
        o_new.pseudorange = o.pseudorange - iono - tropo;
        corrected_obs.push_back(o_new);

        cout << "  G" << o.prn << " el=" << el
            << " iono=" << iono << " tropo=" << tropo << endl;
    }

    PositionSolution sol2 = solveGpsPosition(corrected_obs, false);

    // ─────────────────────────────────────────────────────────────
    // Compare PASS 1 vs PASS 2 (residual metrics)
    // ─────────────────────────────────────────────────────────────
    cout << endl << "==============================" << endl;
    cout << "    BEFORE vs AFTER CORRECTIONS   " << endl;
    cout << "==============================" << endl;
    cout << fixed << setprecision(2);
    cout << "PASS 1 (no atmo):   sigma-0 = " << sol1.sigma0
        << " m, PDOP = " << sol1.dop.PDOP << endl;
    cout << "PASS 2 (with atmo): sigma-0 = " << sol2.sigma0
        << " m, PDOP = " << sol2.dop.PDOP << endl;
    cout << "Sigma-0 reduction: " << (sol1.sigma0 - sol2.sigma0)
        << " m (" << (100.0 * (sol1.sigma0 - sol2.sigma0) / sol1.sigma0)
        << "%)" << endl;

    // ─────────────────────────────────────────────────────────────
    // Print final position in LLA
    // ─────────────────────────────────────────────────────────────
    double final_lat, final_lon, final_h;
    ecefToLla(sol2.x, sol2.y, sol2.z, final_lat, final_lon, final_h);
    cout << endl << "=== FINAL POSITION ===" << endl;
    cout << "  Latitude:  " << final_lat << " deg" << endl;
    cout << "  Longitude: " << final_lon << " deg" << endl;
    cout << "  Height:    " << final_h << " m" << endl;

    // ─────────────────────────────────────────────────────────────
    // Compare against true BILL00USA position (from IGS site log)
    // ─────────────────────────────────────────────────────────────

    double err1 = sqrt(pow(sol1.x - TRUE_X, 2)
        + pow(sol1.y - TRUE_Y, 2)
        + pow(sol1.z - TRUE_Z, 2));

    double err2 = sqrt(pow(sol2.x - TRUE_X, 2)
        + pow(sol2.y - TRUE_Y, 2)
        + pow(sol2.z - TRUE_Z, 2));

    cout << endl << "=== ACCURACY vs TRUE POSITION ===" << endl;
    cout << "PASS 1 (no atmo):   3D error = " << err1 << " m" << endl;
    cout << "PASS 2 (with atmo): 3D error = " << err2 << " m" << endl;
    cout << "Position improvement: " << (err1 - err2) << " m" << endl;

    return 0;
}