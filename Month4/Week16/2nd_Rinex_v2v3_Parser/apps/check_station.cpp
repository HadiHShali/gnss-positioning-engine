// ============================================================
// check_station.cpp -- one-station health check, run BEFORE any differencing
//
// Usage:   check_station <obs_file> [nav_file]
// Example: ./build/Release/check_station.exe data/bill0070.24o data/bill0070.24n
//
// OBS part (needs only the obs file):
//   which reader ran (v2/v3), which codes were picked, approx XYZ,
//   sampling interval, first/last epoch, satellites, dual-frequency share
//
// NAV part (only if a nav file is given):
//   number of ephemerides, unhealthy satellites, and a PARSING sanity check:
//   if a column offset were wrong (e.g. margin 3 vs 4), the orbit values would
//   be wrong -- so we check they look like GPS orbits:
//     semi-major axis A = sqrt_a^2  about 26,560 km   (GPS orbit radius)
//     eccentricity e                below ~0.03       (GPS orbits are near-circular)
//     inclination i0                about 55 deg      (+/- a few degrees)
//     node rate Omega_dot           NEGATIVE, ~ -8e-9 rad/s for every GPS satellite
//     GPS week                      same as the obs file's week
//   and finally: which observed satellites have NO usable (healthy) ephemeris.
//
//   Why the Omega_dot sign test matters: stod() reads the longest valid number
//   at the start of a field and ignores the rest. Reading one column too far
//   RIGHT therefore still gives the right digits -- but drops the minus sign!
//   A, e and i0 are always positive, so only a quantity that is ALWAYS
//   negative (Omega_dot: GPS orbit planes always drift westward) catches it.
//   Reading one column too far LEFT cuts the exponent ("E+03" -> "E+0") and
//   the A check catches that.
//
// Uses only parseRinexAllEpochs (rinex_obs.h) and parseRinexNav /
// findBestEphemeris (satellites.h).
// ============================================================

#include <iostream>      // cout (normal output), cerr (error output), endl
#include <iomanip>       // output formatting: setw, setfill, setprecision, fixed, scientific
#include <set>           // std::set: a sorted collection where every value appears only once
#include <cmath>         // M_PI (needs _USE_MATH_DEFINES on MSVC -- set in CMakeLists.txt)
#include <algorithm>     // std::min, std::max
#include "rinex_obs.h"   // parseRinexAllEpochs, ObsEpoch, ObsRecord, RinexObsHeader
#include "satellites.h"  // parseRinexNav, findBestEphemeris, GpsEphemeris

using namespace std;     // lets us write cout instead of std::cout, set instead of std::set, ...

// ------------------------------------------------------------
// printTime: print one epoch's date, time and GPS time on one line, e.g.
//   2024-01-07 00:00:00.0  (GPS week 2296, sow 0.0)
// 'static' = only visible inside this file.
// 'const ObsEpoch& e' = read-only reference: no copy is made, nothing is changed.
// ------------------------------------------------------------
static void printTime(const ObsEpoch& e) {
    cout << e.year << "-"                         // year, then a dash
        << setfill('0') << setw(2) << e.month    // setfill('0') + setw(2): pad to 2 digits with zeros -> "01"
        << "-" << setw(2) << e.day               // setw applies to the NEXT value only, so repeat it
        << " " << setw(2) << e.hour              // hour, 2 digits
        << ":" << setw(2) << e.minute            // minute, 2 digits
        << ":" << setw(4) << fixed << setprecision(1) << e.second   // seconds: 4 wide, 1 decimal -> "00.0"
        << setfill(' ')                          // switch the fill character back to blanks
        << "  (GPS week " << e.gps_week          // GPS week computed by the parser (2296)
        << ", sow " << e.t_gps << ")";           // seconds of week (0 at Sunday 00:00)
}

// ------------------------------------------------------------
// main receives everything typed on the command line:
//   argc = argument count:  how many words, INCLUDING the program name
//   argv = argument vector: an array of C-strings (char*), one per word
//
// Example:
//   ./build/Release/check_station.exe data/lchs...MO.rnx data/lchs...MN.rnx
//   argc    = 3
//   argv[0] = "./build/Release/check_station.exe"   (the program itself)
//   argv[1] = "data/lchs...MO.rnx"                  (obs file, required)
//   argv[2] = "data/lchs...MN.rnx"                  (nav file, optional)
//   argv[3] = nullptr                               (marks the end of the list)
//
// The shell splits the line at spaces, so a path that contains spaces
// must be put in quotes to arrive as ONE argument.
// ------------------------------------------------------------
int main(int argc, char* argv[]) {      // argc: argument count, argv: argument vector
    if (argc < 2) {                     // only the program name was typed: no obs file given
        // Reading argv[1] now would read past the end of the array (undefined
        // behavior), so print how to use the program and stop instead.
        // argv[0] is used so the message always shows the program's real name.
        cerr << "Usage: " << argv[0] << " <obs_file> [nav_file]" << endl;
        // <...> marks a required argument, [...] an optional one (common convention)
        return 1;                       // non-zero exit code = "something went wrong" (0 = success)
    }

    // ==========================================================
    // OBS part
    // ==========================================================
    cout << "=== STATION CHECK: " << argv[1] << " ===" << endl;    // title line with the obs file name

    RinexObsHeader h;                                   // empty header summary; the parser fills it in
    const auto epochs = parseRinexAllEpochs(argv[1], &h);   // read the whole obs file (v2 or v3).
    //   &h = address of h -> the header comes back in h.
    //   'auto' = the compiler works out the type:
    //   vector<ObsEpoch>. 'const' = we won't change it.
    //   The parser itself prints the version line and
    //   the "Resolved indices" line.
    if (epochs.empty()) {                               // nothing parsed: missing file, compressed, CRINEX, ...
        cerr << "No epochs parsed -- see the messages above." << endl;   // the parser already said why
        return 1;                                       // stop with an error code
    }

    // ---- Statistics over all epochs ----
    set<int> prns;                                      // every GPS PRN seen at least once (no duplicates, sorted)
    size_t n_rec = 0;                                   // total satellite records over all epochs
    //   (size_t = unsigned integer type used for sizes/counts)
    size_t n_dual = 0;                                  // ... of which have P1, P2, L1 and L2 all present
    size_t max_sats = 0;                                // the largest number of satellites in one epoch

    for (const auto& ep : epochs) {                     // loop over every epoch (read-only reference, no copy)
        max_sats = max(max_sats, ep.gps_records.size());    // keep the larger of: best so far, this epoch's count
        for (const auto& r : ep.gps_records) {          // loop over every satellite record in this epoch
            prns.insert(r.prn);                         // remember this PRN (set ignores it if already there)
            n_rec++;                                    // count the record
            if (r.has_c2 && r.has_phase_l1 && r.has_phase_l2)   // P1 is always present (parser rule), so
                n_dual++;                               //   these three complete the dual-frequency set
        }
    }

    // Real sampling = time between the first two epochs. The header's INTERVAL
    // line is optional and can be missing or wrong, so we also measure it.
    const double dt = (epochs.size() > 1)               // condition ? value-if-true : value-if-false
        ? epochs[1].t_gps - epochs[0].t_gps //   two or more epochs: their time difference
        : 0.0;                              //   only one epoch: no interval to measure

    cout << "GPS obs types   : ";                       // label
    for (const auto& t : h.gps_obs_types) cout << t << " ";   // every GPS code from the header, in file order
    cout << endl;                                       // end the line

    cout << fixed << setprecision(3);                   // fixed-point notation, 3 decimals (mm) from here on
    cout << "Approx XYZ (m)  : " << h.approx_x << "  "  // header position X
        << h.approx_y << "  " << h.approx_z << endl;   // ... Y and Z (ECEF, meters)

    cout << setprecision(1);                            // 1 decimal from here on
    cout << "Interval        : header " << h.interval   // what the header claims
        << " s | data " << dt << " s" << endl;         // what the data actually shows

    cout << "First epoch     : "; printTime(epochs.front()); cout << endl;   // .front() = first element
    cout << "Last epoch      : "; printTime(epochs.back());  cout << endl;   // .back()  = last element
    cout << "Epochs          : " << epochs.size() << endl;                   // 5760 = full day at 15 s

    cout << "GPS PRNs seen   : " << prns.size()         // number of distinct satellites
        << "  (max " << max_sats << " per epoch, mean "
        << double(n_rec) / epochs.size() << ")" << endl;   // double(...) avoids integer division

    cout << "Dual-freq recs  : "
        << (n_rec ? 100.0 * n_dual / n_rec : 0.0)      // percentage; the ?: guards against dividing by 0
        << " %  (P1 + P2 + L1 + L2 all present)" << endl;

    if (argc < 3) return 0;                             // no nav file was given: obs check done, success

    // ==========================================================
    // NAV part
    // ==========================================================
    cout << endl << "--- NAV: " << argv[2] << " ---" << endl;   // blank line, then a section title

    const auto ephs = parseRinexNav(argv[2]);           // read every GPS ephemeris (v2 or v3);
    //   prints the "NAV file: RINEX ..." line itself
    if (ephs.empty()) {                                 // nothing parsed
        cerr << "No ephemerides parsed -- see the messages above." << endl;
        return 1;
    }

    // ---- One pass over all records: counts, unhealthy PRNs, value ranges ----
    set<int> nav_prns, unhealthy;                       // PRNs that have records / PRNs with sv_health != 0
    double a_min = 1e18, a_max = 0, e_max = 0;          // semi-major axis range [km], largest eccentricity
    //   (min starts huge and max starts small, so the
    //    first record always replaces them)
    double i_min = 1e18, i_max = 0;                     // inclination range [deg]
    double od_min = 1e18, od_max = -1e18;               // node rate Omega_dot range [rad/s] (values are negative,
    //   so the max starts at a hugely negative number)
    set<int> weeks;                                     // every GPS week found in the records

    for (const auto& e : ephs) {                        // loop over every ephemeris record
        nav_prns.insert(e.prn);                         // remember the PRN
        if (e.sv_health != 0.0) unhealthy.insert(e.prn);   // health word not 0 -> satellite flagged unusable
        const double a_km = e.sqrt_a * e.sqrt_a / 1000.0;   // A = (sqrt A)^2, meters -> km
        const double i_deg = e.i0 * 180.0 / M_PI;           // inclination: radians -> degrees
        a_min = min(a_min, a_km);  a_max = max(a_max, a_km);    // widen the A range
        e_max = max(e_max, e.e);                                // keep the largest eccentricity
        i_min = min(i_min, i_deg); i_max = max(i_max, i_deg);   // widen the inclination range
        od_min = min(od_min, e.Omega_dot);                      // widen the Omega_dot range
        od_max = max(od_max, e.Omega_dot);
        weeks.insert(static_cast<int>(e.gps_week));     // gps_week is stored as a double -> convert to int
    }

    cout << "Ephemerides     : " << ephs.size()         // total records
        << " for " << nav_prns.size() << " PRNs" << endl;   // ... and distinct satellites

    cout << "Unhealthy PRNs  : ";
    if (unhealthy.empty()) cout << "none";              // no flagged satellites
    for (int p : unhealthy)                             // otherwise list them as G01, G13, ...
        cout << "G" << setw(2) << setfill('0') << p << setfill(' ') << " ";
    cout << endl;

    // ---- Parsing sanity: would a wrong column offset have slipped through? ----
    cout << setprecision(1);                            // 1 decimal for km and degrees
    cout << "Semi-major axis : " << a_min << " - " << a_max << " km   (GPS: ~26,560)  "
        << ((a_min > 26000 && a_max < 27100) ? "OK" : "<-- CHECK PARSING") << endl;   // GPS orbit size

    cout << setprecision(4);                            // 4 decimals for eccentricity
    cout << "Eccentricity max: " << e_max << "            (GPS: < ~0.03)  "
        << ((e_max < 0.03) ? "OK" : "<-- CHECK PARSING") << endl;                      // near-circular

    cout << setprecision(1);
    cout << "Inclination     : " << i_min << " - " << i_max << " deg      (GPS: ~55)  "
        << ((i_min > 50 && i_max < 60) ? "OK" : "<-- CHECK PARSING") << endl;          // ~55 deg planes

    cout << scientific << setprecision(2);              // scientific notation for tiny numbers: -8.16e-09
    cout << "Node rate       : " << od_min << " - " << od_max << " rad/s  (GPS: all < 0)  "
        << ((od_max < 0) ? "OK" : "<-- CHECK PARSING (lost minus signs?)") << endl;    // must all be negative

    cout << fixed << setprecision(1);                   // back to fixed-point notation
    cout << "GPS week(s)     : ";
    for (int w : weeks) cout << w << " ";               // usually one week (2296); two if records cross Sunday
    cout << "  (obs file: " << epochs.front().gps_week << ")  "
        << (weeks.count(epochs.front().gps_week)       // .count() = 1 if the set contains it, else 0
            ? "OK" : "<-- MISMATCH: different week?") << endl;

    // ---- Coverage: does every observation have a usable ephemeris? ----
    // findBestEphemeris returns nullptr if a satellite has no healthy record
    // within 2 h of Toe. Checking every record (not one mid-file time) matters:
    // a satellite seen only in the morning has no record near noon, and that
    // would be a false alarm.
    size_t n_no_eph = 0;                                // observation records without a usable ephemeris
    set<int> prns_no_eph;                               // ... and which satellites they belong to
    for (const auto& ep : epochs)                       // every epoch
        for (const auto& r : ep.gps_records)            // every satellite record in it
            if (!findBestEphemeris(ephs, r.prn, ep.t_gps)) {   // nullptr: none, unhealthy, or > 2 h old
                n_no_eph++;
                prns_no_eph.insert(r.prn);
            }
    cout << "No usable eph   : " << n_no_eph << " of " << n_rec << " records";
    if (!prns_no_eph.empty()) {                         // list the satellites affected
        cout << "  (";
        for (int p : prns_no_eph) cout << "G" << setw(2) << setfill('0') << p << setfill(' ') << " ";
        cout << ")";
    }
    cout << endl;

    return 0;                                           // everything ran: exit code 0 = success
}