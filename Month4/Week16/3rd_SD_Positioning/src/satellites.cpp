#include "satellites.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>

using namespace std;

const double MU_EARTH = 3.986005e14;
const double OMEGA_EARTH = 7.2921151467e-5;

// Helper functions (static = private to this file)
static double solveKepler(double M, double e) {
    double E = M;
    for (int i = 0; i < 20; i++) {
        double E_new = M + e * sin(E);
        if (abs(E_new - E) < 1e-12) return E_new;
        E = E_new;
    }
    return E;
}

static string fortranToCpp(const string& s) {
    string out = s;
    for (char& c : out) if (c == 'D' || c == 'd') c = 'E';
    return out;
}

static double parseField(const string& line, size_t start, size_t len = 19) {
    if (start >= line.length()) return 0.0;
    string field = fortranToCpp(line.substr(start, len));
    try { return stod(field); }
    catch (...) { return 0.0; }
}

// ------------------------------------------------------------
// readOrbitLine: read the 4 numbers of one BROADCAST ORBIT line.
//
// Each orbit line holds 4 values, each 19 characters wide (Fortran D19.12),
// written with a 'D' exponent, e.g. " 6.400000000000D+01". parseField()
// turns the 'D' into 'E' and converts the text to a double.
//
// The ONLY layout difference between versions is the left margin:
//   RINEX 2:  3X, 4D19.12  -> 3 blanks; values start at 3, 22, 41, 60
//   RINEX 3:  4X, 4D19.12  -> 4 blanks; values start at 4, 23, 42, 61
// So the caller passes the margin as 'off' (3 or 4), and the values sit
// at off, off+19, off+38, off+57.
//
// Example, bill0070.24n line 10 (RINEX 2, PRN 1, BROADCAST ORBIT 5):
//   "   -2.542963067476D-10 1.000000000000D+00 2.296000000000D+03 0.000000000000D+00"
//    ^^^ 3 blanks -> IDOT at 3, L2 codes at 22, GPS week at 41, L2P flag at 60
//
// The last orbit line often has only 2 values; parseField() returns 0.0
// for a field that starts past the end of the line, so that is safe.
// ------------------------------------------------------------
static void readOrbitLine(const string& line, int off,
	double& v1, double& v2, double& v3, double& v4) {
	v1 = parseField(line, off);         // value 1 at the margin
	v2 = parseField(line, off + 19);    // each value is 19 characters wide
	v3 = parseField(line, off + 38);
	v4 = parseField(line, off + 57);    // may be absent on the last line -> 0.0
}

// ======================================================================
// parseRinexNav: read every GPS broadcast ephemeris from a RINEX
// navigation file, version 2.xx or 3.xx.
//
// What both versions share: after END OF HEADER, each satellite record is
// 8 lines -- 1 first line (satellite, Toc epoch, 3 clock terms) + 7
// BROADCAST ORBIT lines (4 values each, the last one usually 2). The 29
// values and their order are IDENTICAL in both versions; only the layout
// of the first line and the orbit-line margin differ.
//
// What only RINEX 3 has: several systems in one file ('M' = mixed).
// Non-GPS records must be skipped with the RIGHT length:
//   GLONASS (R), SBAS (S):                  first line + 3 orbit lines
//   Galileo, BeiDou, QZSS, NavIC (E,C,J,I): first line + 7 orbit lines
// The old parser always skipped 7 lines after a non-GPS first line. With
// GLONASS records in between it lost track: in a test file with GLONASS
// records interleaved it kept only 81 of 162 GPS records. It worked on the
// BRDC file only because of the order of the records in that file.
// ======================================================================
vector<GpsEphemeris> parseRinexNav(const string& filename) {
    vector<GpsEphemeris> ephemerides;                   // result; stays empty on any error
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "ERROR: Cannot open " << filename << endl;
        return ephemerides;
    }

    // ---- Line 1: RINEX VERSION / TYPE ----
    //   index 0-8   version (F9.2)      "     2.11" / "     3.04"
    //   index 20    file type           'N' = navigation
    //   index 40    satellite system    blank = GPS (v2), 'G' = GPS, 'M' = mixed (v3)
    string line;
    if (!getline(file, line)) {
        cerr << "ERROR: " << filename << " is empty" << endl;
        return ephemerides;
    }
    if (!line.empty() && line[0] == '\x1f') {            // .gz and .Z files start with byte 0x1F
        cerr << "ERROR: " << filename << " is still compressed (.gz/.Z) -- run: gzip -d <file>" << endl;
        return ephemerides;
    }
    if (line.find("RINEX VERSION / TYPE") == string::npos) {
        cerr << "ERROR: " << filename << ": line 1 is not 'RINEX VERSION / TYPE'" << endl;
        return ephemerides;
    }
    const double version = parseField(line, 0, 9);      // 2.11 or 3.04
    const int major = static_cast<int>(version);        // 2 or 3
    const char ftype = (line.size() > 20) ? line[20] : ' ';
    const char fsys = (line.size() > 40) ? line[40] : ' ';

    if (ftype != 'N') {                                  // e.g. 'O': an obs file passed by mistake
        cerr << "ERROR: not a GPS/mixed nav file (type '" << ftype << "')" << endl;
        return ephemerides;
    }
    if (major == 3 && fsys != 'G' && fsys != 'M') {      // e.g. 'E': a Galileo-only nav file
        cerr << "ERROR: RINEX 3 nav file for system '" << fsys << "', no GPS records" << endl;
        return ephemerides;
    }
    if (major < 2 || major > 3) {                        // RINEX 4 nav uses "> EPH" record headers
        cerr << "ERROR: NAV RINEX " << version << " not supported" << endl;
        return ephemerides;
    }
    string ver = line.substr(0, 9);                      // print the version exactly as written
    ver.erase(0, ver.find_first_not_of(' '));            //   "     2.11" -> "2.11"
    cout << "NAV file: RINEX " << ver << " -> v" << major << " reader" << endl;

    // ---- Skip the rest of the header ----
    // (IONOSPHERIC CORR / ION ALPHA, TIME SYSTEM CORR, LEAP SECONDS, ... are
    //  not needed for orbits; parseKlobucharFromNav reads the iono lines.)
    bool header_done = false;
    while (getline(file, line)) {
        if (line.find("END OF HEADER") != string::npos) { header_done = true; break; }
    }
    if (!header_done) {
        cerr << "ERROR: No END OF HEADER found" << endl;
        return ephemerides;
    }

    const int off = (major == 2) ? 3 : 4;                // orbit-line margin (see readOrbitLine)

    // ---- Records ----
    while (getline(file, line)) {
        if (line.find_first_not_of(" \r") == string::npos) continue;   // blank line

        GpsEphemeris eph = {};                           // = {} sets every field to 0

        if (major == 2) {
            // RINEX 2 first line (file is GPS-only, so there is no system letter):
            //   " 1 24  1  7  0  0  0.0 1.654480583966D-04 9.094947017729D-13 0.000000000000D+00"
            //   index 0-1  PRN          (I2)      " 1"
            //         3-4  year         (I2.2)    "24"  -> 2024 (80-99 -> 19xx)
            //         6-7  month        (I2)      " 1"
            //         9-10 day          (I2)      " 7"
            //        12-13 hour         (I2)      " 0"
            //        15-16 minute       (I2)      " 0"
            //        17-21 second       (F5.1)    "  0.0"
            //        22 / 41 / 60  clock bias, drift, drift rate (3 x D19.12)
            eph.prn = static_cast<int>(parseField(line, 0, 2));
            const int yy = static_cast<int>(parseField(line, 3, 2));
            eph.year = (yy < 80) ? 2000 + yy : 1900 + yy;
            eph.month = static_cast<int>(parseField(line, 6, 2));
            eph.day = static_cast<int>(parseField(line, 9, 2));
            eph.hour = static_cast<int>(parseField(line, 12, 2));
            eph.minute = static_cast<int>(parseField(line, 15, 2));
            eph.second = parseField(line, 17, 5);
            eph.clk_bias = parseField(line, 22);
            eph.clk_drift = parseField(line, 41);
            eph.clk_drift_rate = parseField(line, 60);
        }
        else {
            // RINEX 3 first line:
            //   "G01 2024 01 07 00 00 00 1.654480583966E-04 9.094947017729E-13 0.000000000000E+00"
            //   index 0    system letter  (A1)      'G'
            //         1-2  PRN            (I2.2)    "01"
            //         4-7  year           (I4)      "2024"
            //         9-10 month, 12-13 day, 15-16 hour, 18-19 minute, 21-22 second (I2 each)
            //        23 / 42 / 61  clock bias, drift, drift rate (3 x D19.12)
            const char sys = line[0];
            if (sys != 'G') {
                // Not GPS: skip the rest of this record with the correct length.
                const int skip = (sys == 'R' || sys == 'S') ? 3 : 7;
                for (int i = 0; i < skip; i++) if (!getline(file, line)) break;
                continue;
            }
            eph.prn = static_cast<int>(parseField(line, 1, 2));
            eph.year = static_cast<int>(parseField(line, 4, 4));
            eph.month = static_cast<int>(parseField(line, 9, 2));
            eph.day = static_cast<int>(parseField(line, 12, 2));
            eph.hour = static_cast<int>(parseField(line, 15, 2));
            eph.minute = static_cast<int>(parseField(line, 18, 2));
            eph.second = parseField(line, 21, 2);
            eph.clk_bias = parseField(line, 23);
            eph.clk_drift = parseField(line, 42);
            eph.clk_drift_rate = parseField(line, 61);
        }

        // The 7 BROADCAST ORBIT lines -- same content and order in v2 and v3.
        // This is the GPS LNAV broadcast ephemeris (IS-GPS-200): the orbit and
        // clock model the Control Segment uploads to each satellite.
        // Angles are in RADIANS (the satellite broadcasts semicircles; RINEX
        // converts). Example values: G01, week 2296, Toe = Sunday 00:00.
        //
        // (The first line, read above, holds the CLOCK model:
        //    Toc    reference time of the clock model         2024-01-07 00:00:00
        //    af0    clock bias at Toc [s]                     1.6545e-4 s = 165 us ~ 49.6 km of range!
        //    af1    clock drift [s/s]                         9.09e-13
        //    af2    clock drift rate [s/s^2]                  0
        //  dt_sv = af0 + af1*(t-Toc) + af2*(t-Toc)^2 + relativistic term.)
        //
        // How the orbit parameters fit together (what computeSatPosECEF does):
        //   * 6 Kepler elements describe an ideal ellipse at Toe:
        //       sqrt(A), e, i0, OMEGA0, omega, M0
        //   * 3 rates move it forward in time:
        //       delta_n (mean motion), OMEGA_DOT (node), IDOT (inclination)
        //   * 6 harmonic corrections add the real perturbations, each as
        //     sin/cos of TWICE the argument of latitude (Phi = true anomaly + omega):
        //       radius        dr = Crs*sin(2Phi) + Crc*cos(2Phi)   [m]
        //       along-track   du = Cus*sin(2Phi) + Cuc*cos(2Phi)   [rad]
        //       inclination   di = Cis*sin(2Phi) + Cic*cos(2Phi)   [rad]
        // ==================================================================

        // ---- Orbit 1 ----------------------------------------------------
        //   IODE     Issue Of Data, Ephemeris: counter that changes with
        //            every new orbit upload                          64
        //   Crs      sine correction to the orbit RADIUS [m]         -68.06 m
        //   delta_n  mean-motion correction [rad/s]:
        //            n = sqrt(mu / A^3) + delta_n                    4.01e-9 (n ~ 1.458e-4)
        //   M0       MEAN ANOMALY at Toe: where the satellite is
        //            along its orbit at the reference time [rad]     0.404 rad = 23.1 deg
        if (!getline(file, line)) break;
        readOrbitLine(line, off, eph.IODE, eph.Crs, eph.delta_n, eph.M0);

        // ---- Orbit 2 ----------------------------------------------------
        //   Cuc      cosine correction to the ARGUMENT OF LATITUDE
        //            (along-track) [rad]                             -3.45e-6 rad ~ -92 m
        //   e        ECCENTRICITY (0 = circle)                       0.0131
        //   Cus      sine correction to the argument of latitude
        //            [rad]                                           4.36e-6 rad ~ 116 m
        //   sqrt_a   square root of the SEMI-MAJOR AXIS [sqrt(m)]    5154.02 -> A = 26,564 km
        //            (with e = 0.0131: r between ~26,216 and ~26,912 km)
        if (!getline(file, line)) break;
        readOrbitLine(line, off, eph.Cuc, eph.e, eph.Cus, eph.sqrt_a);

        // ---- Orbit 3 ----------------------------------------------------
        //   toe      reference time of the EPHEMERIS
        //            [s of GPS week]                                 0 s = Sunday 00:00
        //   Cic      cosine correction to the INCLINATION [rad]      -1.27e-7 rad ~ -3.4 m
        //   Omega0   LONGITUDE OF THE ASCENDING NODE at the START
        //            OF THE GPS WEEK -- not at Toe! [rad]            -1.734 rad = -99.4 deg
        //            (this is why computeSatPosECEF subtracts OMEGA_EARTH * toe)
        //   Cis      sine correction to the inclination [rad]        1.01e-7 rad ~ 2.7 m
        if (!getline(file, line)) break;
        readOrbitLine(line, off, eph.toe, eph.Cic, eph.Omega0, eph.Cis);

        // ---- Orbit 4 ----------------------------------------------------
        //   i0        INCLINATION of the orbit plane at Toe [rad]    0.990 rad = 56.7 deg
        //   Crc       cosine correction to the orbit RADIUS [m]      310.7 m
        //   omega     ARGUMENT OF PERIGEE: where in the plane the
        //             closest point to Earth lies [rad]              0.999 rad = 57.2 deg
        //   Omega_dot drift RATE of the ascending node [rad/s]       -8.16e-9 = -0.040 deg/day
        if (!getline(file, line)) break;
        readOrbitLine(line, off, eph.i0, eph.Crc, eph.omega, eph.Omega_dot);

        // ---- Orbit 5 ----------------------------------------------------
        //   i_dot     RATE of the inclination [rad/s]                -2.54e-10
        //   L2_codes  code on L2: 1 = P code, 2 = C/A code           1
        //   gps_week  GPS week that Toe belongs to; RINEX writes a
        //             continuous count, NOT modulo 1024              2296
        //   L2P_flag  0 = navigation data is on the L2 P code        0
        if (!getline(file, line)) break;
        readOrbitLine(line, off, eph.i_dot, eph.L2_codes, eph.gps_week, eph.L2P_flag);

        // ---- Orbit 6 ----------------------------------------------------
        //   sv_accuracy  URA: predicted range accuracy [m]           2.8 m
        //   sv_health    0 = healthy, anything else = DO NOT USE    63 -> G01 is UNHEALTHY
        //                (G01 is also missing from the BILL, NMKM and LCHS obs tables.
        //                 TODO: skip records with sv_health != 0 in findBestEphemeris)
        //   TGD          group delay L1 vs L2 P(Y) [s]               5.12e-9 s ~ 1.54 m
        //                single-frequency L1 users must apply it; the broadcast
        //                clock refers to the ionosphere-free combination, so the
        //                dual-frequency IF solution (Week 13) correctly ignores it
        //   IODC         Issue Of Data, Clock: must match IODE for a
        //                consistent orbit + clock set                64 (= IODE)
        if (!getline(file, line)) break;
        readOrbitLine(line, off, eph.sv_accuracy, eph.sv_health, eph.TGD, eph.IODC);

        // ---- Orbit 7 ----------------------------------------------------
        //   trans_time    when the satellite SENT this message
        //                 [s of GPS week]                            -7182 s
        //                 negative = before the week began:
        //                 604800 - 7182 = Saturday 22:00:18, ~2 h before Toe
        //   fit_interval  how long the orbit fit is valid [hours]   4 h -> Toe +/- 2 h
        //   (3rd and 4th slots are spares, usually absent -> parseField returns 0)
        if (!getline(file, line)) break;
        double spare1, spare2;
        readOrbitLine(line, off, eph.trans_time, eph.fit_interval, spare1, spare2);

        ephemerides.push_back(eph);                      // one complete GPS record
    }
    return ephemerides;
}

// - computeSatPosECEF()
SatPosition computeSatPosECEF(const GpsEphemeris& eph, double t_gps)
{
	SatPosition pos;
	double a = eph.sqrt_a * eph.sqrt_a;
	double n = sqrt(MU_EARTH / (a * a * a)) + eph.delta_n;
	double tk = t_gps - eph.toe;
	if (tk > 302400.0) tk -= 604800.0;
	if (tk < -302400.0) tk += 604800.0;
	double Mk = eph.M0 + n * tk;
	double Ek = solveKepler(Mk, eph.e);
	double nuk = atan2(sqrt(1 - eph.e * eph.e) * sin(Ek), cos(Ek) - eph.e);
	double phi_k = nuk + eph.omega;
	double sin_2phi = sin(2 * phi_k), cos_2phi = cos(2 * phi_k);
	double u_k = phi_k + eph.Cuc * cos_2phi + eph.Cus * sin_2phi;
	double r_k = a * (1 - eph.e * cos(Ek)) + eph.Crc * cos_2phi + eph.Crs * sin_2phi;
	double i_k = eph.i0 + eph.i_dot * tk + eph.Cic * cos_2phi + eph.Cis * sin_2phi;
	double x_orb = r_k * cos(u_k);
	double y_orb = r_k * sin(u_k);
	double Omega_k = eph.Omega0 + (eph.Omega_dot - OMEGA_EARTH) * tk - OMEGA_EARTH * eph.toe;
	pos.X = x_orb * cos(Omega_k) - y_orb * cos(i_k) * sin(Omega_k);
	pos.Y = x_orb * sin(Omega_k) + y_orb * cos(i_k) * cos(Omega_k);
	pos.Z = y_orb * sin(i_k);
	return pos;
}
// ------------------------------------------------------------
// ------------------------------------------------------------
// findBestEphemeris: choose the ephemeris to use for one satellite at time t.
//
//   ephs   all GPS records from parseRinexNav (healthy AND unhealthy)
//   prn    satellite number
//   t_gps  time of interest [s of GPS week]
//   returns a pointer to the chosen record, or nullptr if there is none
//
// Rule 1 -- health: skip every record whose sv_health is not 0 (the Control
//   Segment's "do not use"). Example: G01 on 2024-01-07, sv_health = 63.
// Rule 2 -- closest Toe: among the healthy records, take the one whose Toe is
//   nearest to t_gps, because the orbit fit is most accurate near Toe.
// Rule 3 -- age limit: reject the winner if it is more than 2 h from t_gps.
//   A GPS ephemeris is fitted over 4 h (fit interval in orbit line 7), i.e.
//   it is meant for Toe +/- 2 h; outside that the orbit error grows quickly.
//   (RTKLIB uses the same 7200 s limit for GPS.)
//
// nullptr means: no record, all records unhealthy, or all too old. The caller
// must skip the satellite at this epoch.
// ------------------------------------------------------------
const GpsEphemeris* findBestEphemeris(const vector<GpsEphemeris>& ephs, int prn, double t_gps)
{
    const double MAX_AGE = 7200.0;           // [s] Toe +/- 2 h (4-hour fit interval)
    const GpsEphemeris* best = nullptr;      // nullptr means "no pointer yet -- empty". We start with nothing found, then update as we discover matches.
    double best_dt = 1e18;                   // time difference of the best match so far. We want the SMALLEST dt, so we start at a deliberately HUGE value.
    for (const auto& e : ephs) {             // read-only reference to each record in the file
        if (e.prn != prn) continue;          // a different satellite
        if (e.sv_health != 0.0) continue;    // unhealthy -> never use it, however close its Toe is
        double dt = abs(e.toe - t_gps);
        if (dt > 302400.0) dt = 604800.0 - dt;           // week rollover: t near Saturday 24:00, Toe near Sunday 00:00
        if (dt < best_dt) { best_dt = dt; best = &e; }   // '&e' gives the address of e in memory
    }
    return (best_dt <= MAX_AGE) ? best : nullptr;        // NEW: closest record, but only if it is fresh enough
}

// compute elevation angle (in degrees) from receiver to satellite
double computeElevation(double rec_x, double rec_y, double rec_z, double sat_x, double sat_y, double sat_z)
{
	// receiver lat lon (approximate from ECEF)
	double r = sqrt(rec_x*rec_x + rec_y*rec_y + rec_z*rec_z);
	double lat = atan2(rec_z, sqrt(rec_x*rec_x+rec_y*rec_y));
	double lon = atan2(rec_y, rec_x);

	// satellite vector relative to receiver
	double dx = sat_x - rec_x;
	double dy = sat_y - rec_y;
	double dz = sat_z - rec_z;

	// rotate ECEF to ENU
	double sin_lat = sin(lat);
	double cos_lat = cos(lat);
	double sin_lon = sin(lon);
	double cos_lon = cos(lon);

	double east = -sin_lon * dx + cos_lon * dy;
	double north = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz;
	double up    =  cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz;

	// elevation angle
	double horizontal = sqrt(east*east + north*north);
	double elev_rad = atan2(up, horizontal);

	return elev_rad * 180.0 / 3.14159265358979;
}

double computeAzimuth (double rec_x, double rec_y, double rec_z, double sat_x, double sat_y, double sat_z)
{
	// same ENU transform as computeElevation, but return azimuth
	double lat = atan2(rec_z, sqrt(rec_x*rec_x + rec_y*rec_y));
	double lon = atan2(rec_y, rec_x);
	double dx = sat_x - rec_x;
	double dy = sat_y - rec_y;
	double dz = sat_z - rec_z;
	double sin_lat = sin(lat), cos_lat = cos(lat);
	double sin_lon = sin(lon), cos_lon = cos(lon);
	double east  = -sin_lon * dx + cos_lon * dy;
	double north = -sin_lat*cos_lon*dx - sin_lat*sin_lon*dy + cos_lat*dz;
    double az_rad = atan2(east, north);  // atan2 returns -pi..pi
    if (az_rad < 0) az_rad += 2 * 3.14159265358979;
    return az_rad * 180.0 / 3.14159265358979;
}


// Parse Klobuchar paramaters from RINEX Nav file header
// ------------------------------------------------------------
// parseKlobucharFromNav: read the GPS Klobuchar ionosphere coefficients
// (alpha0..3, beta0..3) from the HEADER of a RINEX 2 or RINEX 3 nav file.
// Needed only for single-frequency processing; the dual-frequency
// ionosphere-free solution and double differences do not use them.
//
// The same 8 numbers, two layouts (each value 12 chars wide, D12.4):
//
//   RINEX 2  (labels "ION ALPHA" / "ION BETA", format 2X,4D12.4)
//     "    1.1180D-08  1.4900D-08 -5.9600D-08 -5.9600D-08          ION ALPHA"
//      values at 2, 14, 26, 38
//
//   RINEX 3  (label "IONOSPHERIC CORR", type in cols 0-3, format A4,1X,4D12.4)
//     "GPSA   1.9558E-08  0.0000E+00 -5.9605E-08  1.1921E-07       IONOSPHERIC CORR"
//      values at 5, 17, 29, 41      (GPSA = alpha, GPSB = beta; GAL, QZS, ... ignored)
//
// Both lines are optional: e.g. EarthScope's bill0070.24n has neither, so
// valid stays false. Callers must check kp.valid before using the values.
// Units: alpha_n in s/semicircle^n, beta_n in s/semicircle^n.
// ------------------------------------------------------------
KlobucharParams parseKlobucharFromNav(const std::string& filename)
{
    KlobucharParams kp;                                  // alpha/beta start at 0, valid = false
    ifstream file(filename);
    if (!file.is_open()) return kp;

    string line;
    bool found_alpha = false;
    bool found_beta = false;

    while (getline(file, line))                          // header lines only
    {
        if (line.find("END OF HEADER") != string::npos) break;

        // ---- RINEX 3: "IONOSPHERIC CORR", type GPSA / GPSB in columns 0-3 ----
        // find(..., 60): search only from column 60, where header labels live,
        // so text inside a COMMENT line can never match.
        if (line.find("IONOSPHERIC CORR", 60) != string::npos)
        {
            if (line.substr(0, 4) == "GPSA")
            {
                for (int k = 0; k < 4; k++) kp.alpha[k] = parseField(line, 5 + 12 * k, 12);   // 5, 17, 29, 41
                found_alpha = true;
            }
            else if (line.substr(0, 4) == "GPSB")
            {
                for (int k = 0; k < 4; k++) kp.beta[k] = parseField(line, 5 + 12 * k, 12);
                found_beta = true;
            }
        }
        // ---- RINEX 2: separate labels, values start after 2 blanks ----
        else if (line.find("ION ALPHA", 60) != string::npos)
        {
            for (int k = 0; k < 4; k++) kp.alpha[k] = parseField(line, 2 + 12 * k, 12);       // 2, 14, 26, 38
            found_alpha = true;
        }
        else if (line.find("ION BETA", 60) != string::npos)
        {
            for (int k = 0; k < 4; k++) kp.beta[k] = parseField(line, 2 + 12 * k, 12);
            found_beta = true;
        }
    }
    kp.valid = (found_alpha && found_beta);              // usable only if both sets were found
    return kp;
}

// Klobuchar Ionospheric Delay Model (IS-GPS-200)
// Returns delay in meters for L1 frequency
double klobucharIonoDelay (const KlobucharParams& kp, double rec_lat_deg, double rec_lon_deg, double sat_az_deg, double sat_el_deg, double t_gps_sec)
{
	if (!kp.valid) return 0.0; //safe fallback
	    
	const double C_LIGHT = 299792458.0;
    
    // Convert to semi-circles (1 semi-circle = 180 deg = pi rad)
    double phi_u   = rec_lat_deg / 180.0;
    double lambda_u = rec_lon_deg / 180.0;
    double A       = sat_az_deg * M_PI / 180.0;       // radians
    double E       = sat_el_deg / 180.0;              // semi-circles
    
    // STEP 1: Earth-centered angle psi (semi-circles)
    double psi = 0.0137 / (E + 0.11) - 0.022;
	
	// STEP 2: Geodetic latitude of IPP (Ionospheric Pierce Point)
    double phi_i = phi_u + psi * cos(A);
    if (phi_i >  0.416) phi_i =  0.416;
    if (phi_i < -0.416) phi_i = -0.416;
	
	// STEP 3: Geodetic longitude of IPP
    double lambda_i = lambda_u + psi * sin(A) / cos(phi_i * M_PI);

    // STEP 4: Geomagnetic latitude
    double phi_m = phi_i + 0.064 * cos((lambda_i - 1.617) * M_PI);

    // STEP 5: Local time at IPP (seconds)
    double t = 43200.0 * lambda_i + t_gps_sec;
    while (t >= 86400.0) t -= 86400.0;
    while (t < 0.0)      t += 86400.0;

    // STEP 6: Compute amplitude and period from alpha/beta
    double AMP = kp.alpha[0] + kp.alpha[1]*phi_m
               + kp.alpha[2]*phi_m*phi_m + kp.alpha[3]*phi_m*phi_m*phi_m;
    if (AMP < 0.0) AMP = 0.0;
	double PER = kp.beta[0] + kp.beta[1]*phi_m
               + kp.beta[2]*phi_m*phi_m + kp.beta[3]*phi_m*phi_m*phi_m;
    if (PER < 72000.0) PER = 72000.0;

    // STEP 7: Phase
    double x = 2.0 * M_PI * (t - 50400.0) / PER;
    
    // STEP 8: Slant factor F (1 at zenith, ~3 at horizon)
    double F = 1.0 + 16.0 * pow(0.53 - E, 3);
    
    // STEP 9: Compute delay in SECONDS
    double T_iono_sec;
    if (fabs(x) < 1.57) {
        T_iono_sec = F * (5.0e-9 + AMP * (1.0 - x*x/2.0 + x*x*x*x/24.0));
    } else {
        T_iono_sec = F * 5.0e-9;
    }
    
    // Convert delay from SECONDS to METERS
    return T_iono_sec * C_LIGHT;
}

// Saatamonien tropospheric delay model
// uses standard atmosphere (1013.25 mb, 15 deg c, 50% humidity)
// returns delay in meters at the satellite elevation angle
double saastamoinenTropoDelay(
    double rec_lat_deg, double rec_height_m, double sat_el_deg)
{
    // Safety: invalid elevation -> no correction
    if (sat_el_deg < 1.0) return 0.0;
    
    // Standard atmosphere at sea level
    const double P0 = 1013.25;  // mbar
    const double T0 = 288.15;   // K
    const double RH = 0.5;       // relative humidity
    
    // Adjust for receiver height
    double P = P0 * pow(1.0 - 2.26e-5 * rec_height_m, 5.225);
    double T = T0 - 0.0065 * rec_height_m;
    
    // Saturation vapor pressure (Magnus equation), in mbar
    double T_C = T - 273.15;
    double es  = 6.108 * exp(17.27 * T_C / (T_C + 237.3));
    double e   = RH * es;  // actual vapor pressure
    
    // Convert latitude to radians
    double phi = rec_lat_deg * M_PI / 180.0;
    double h_km = rec_height_m / 1000.0;
    
    // Zenith dry (hydrostatic) delay
    double T_dry_zenith = 0.0022768 * P
                        / (1.0 - 0.00266 * cos(2.0 * phi)
                              - 0.00028 * h_km);
    
    // Zenith wet delay
    double T_wet_zenith = 0.0022768 * (1255.0 / T + 0.05) * e;
    
    // Apply slant factor 1/sin(elev)
    double sin_el = sin(sat_el_deg * M_PI / 180.0);
    double T_slant = (T_dry_zenith + T_wet_zenith) / sin_el;
    
    return T_slant;
}

