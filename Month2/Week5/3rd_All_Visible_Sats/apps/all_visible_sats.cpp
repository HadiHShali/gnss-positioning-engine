// ====================================================//
//				All Visible Satellites                 //
// ====================================================//
//Computes positions of all visible GPS satellites from a receiver location

//SECTION0: Headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <iomanip>
#include <algorithm>
using namespace std;


// Constants
const double MU_EARTH = 3.986005e14;
const double OMEGA_EARTH = 7.2921151467e-5;
const double PI = 3.141592653589793;
const double DEG2RAD = PI / 180.0;
const double RAD2DEG = 180.0 / PI;


// WGS84 ellipsoid constants
const double WGS84_A = 6378137.0;           // semi-major axis (m)
const double WGS84_F = 1.0 / 298.257223563;   // flattening
const double WGS84_E2 = 2 * WGS84_F - WGS84_F * WGS84_F;  // eccentricity squared


// Data Structures
struct GpsEphemeris 
{
    int prn;
    int year, month, day, hour, minute;
    double second;
    double clk_bias, clk_drift, clk_drift_rate;
    double IODE, Crs, delta_n, M0;
    double Cuc, e, Cus, sqrt_a;
    double toe, Cic, Omega0, Cis;
    double i0, Crc, omega, Omega_dot;
    double i_dot, L2_codes, gps_week, L2P_flag;
    double sv_accuracy, sv_health, TGD, IODC;
    double trans_time, fit_interval;
};

struct SatPosition { double X, Y, Z; };

struct VisibleSat 
{
    int prn;
    double X, Y, Z;       // ECEF position (m)
    double elevation_deg; // angle above horizon
    double azimuth_deg;   // compass direction
    double range_m;       // distance from receiver
};

// helper functions to parse the RINEX Nav file
// 1st: change the D letter to E in the files
string fortranToCpp(const string& s)
{
    string out = s;
    for (char& c : out) if (c == 'D' || c == 'd') c = 'E' ;
    return out;
}
// 2nd: extracts a fixed-width numeric field from a text line by giving the begining and the length address
double parseField(const string& line, size_t start, size_t len = 19)
{
    if (start >= line.length()) { return 0.0; };
    string field = fortranToCpp(line.substr(start, len));
    try { return stod(field); }
    catch (...) { return 0.0; }
}

// 3rd: extracts four fixed-width numeric values from a RINEX navigation file line
// Without '&' in 'double& v1', the function would receive a copy of the variable, and any changes made
// inside the function would not affect the original variable.
void readOrbitLine(const string& line, double& v1, double& v2, double& v3, double& v4)
{
    v1 = parseField(line, 4);
    v2 = parseField(line, 23);
    v3 = parseField(line, 42);
    v4 = parseField(line, 61);

}

// 4th: RINEX Nav Parser
// vector<GpsEphemeris> means a dynamic array whose elements are of type GpsEphemeris.
// Since GpsEphemeris is a struct, each element in the vector is one complete GpsEphemeris object.

vector<GpsEphemeris> parseRinexNav(const string& filename)
{
    vector<GpsEphemeris> ephs;
    ifstream file(filename);
    if (!file.is_open()) { cerr << "Cannot open file" << endl; return ephs; }

    string line;
    bool header_done = false;
    while (getline(file, line))
    {
        if (line.find("END OF HEADER") != string::npos)
        {
            header_done = true; break;
        }
    }


    if (!header_done) return ephs;

    while (getline(file, line))
    {
        if (line.length() < 23) continue;
        if (line[0] != 'G')
        {
            for (int i = 0; i < 7; i++) if (!getline(file, line)) break;
            continue;
        }



        GpsEphemeris e = {};

        // from here ... 
        try {
            e.prn = stoi(line.substr(1, 2));
            e.year = stoi(line.substr(4, 4));
            e.month = stoi(line.substr(9, 2));
            e.day = stoi(line.substr(12, 2));
            e.hour = stoi(line.substr(15, 2));
            e.minute = stoi(line.substr(18, 2));
            e.second = parseField(line, 21, 2);
        }
        catch (...) { continue; }
        e.clk_bias = parseField(line, 23);
        e.clk_drift = parseField(line, 42);
        e.clk_drift_rate = parseField(line, 61);
        // ... to here is parsing the Line 0. 

        // from here ...
        if (!getline(file, line)) break;
        readOrbitLine(line, e.IODE, e.Crs, e.delta_n, e.M0);
        if (!getline(file, line)) break;
        readOrbitLine(line, e.Cuc, e.e, e.Cus, e.sqrt_a);
        if (!getline(file, line)) break;
        readOrbitLine(line, e.toe, e.Cic, e.Omega0, e.Cis);
        if (!getline(file, line)) break;
        readOrbitLine(line, e.i0, e.Crc, e.omega, e.Omega_dot);
        if (!getline(file, line)) break;
        readOrbitLine(line, e.i_dot, e.L2_codes, e.gps_week, e.L2P_flag);
        if (!getline(file, line)) break;
        readOrbitLine(line, e.sv_accuracy, e.sv_health, e.TGD, e.IODC);
        if (!getline(file, line)) break;
        double spare1, spare2;
        readOrbitLine(line, e.trans_time, e.fit_interval, spare1, spare2);

        // ... to here is parsing the Line 1-7
        ephs.push_back(e);
    }
    return ephs;
}


// ── SATELLITE POSITION (from Day 2) ──────────────────────────────────────
double solveKepler(double M, double e) {
    double E = M;
    for (int i = 0; i < 20; i++) {
        double E_new = M + e * sin(E);
        if (abs(E_new - E) < 1e-12) return E_new;
        E = E_new;
    }
    return E;
}


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


// New Today: LLA(lat lon alt) to ECEF
SatPosition llaToEcef(double lat_deg, double lon_deg, double height_m)
{
    double phi = lat_deg * DEG2RAD;
    double lam = lon_deg * DEG2RAD;
    double sin_phi = sin(phi), cos_phi = cos(phi);
    double sin_lam = sin(lam), cos_lam = cos(lam);
    double N = WGS84_A / sqrt(1 - WGS84_E2 * sin_phi * sin_phi);
    SatPosition pos;
    pos.X = (N + height_m) * cos_phi * cos_lam;
    pos.Y = (N + height_m) * cos_phi * sin_lam;
    pos.Z = (N * (1 - WGS84_E2) + height_m) * sin_phi;
    return pos;
}

// New today: Compute Azimuth, and Elevation angle
void computeElevAzim(double sat_x, double sat_y, double sat_z,
                    double rec_x, double rec_y, double rec_z,
                    double rec_lat_deg, double rec_lon_deg,
                    double& el_deg, double& az_deg, double& range_m)
{
    // Relative vector (satellite - receiver, in ECEF)
    double dx = sat_x - rec_x;
    double dy = sat_y - rec_y;
    double dz = sat_z - rec_z;


    // Convert receiver lat/lon to radians
    double phi = rec_lat_deg * DEG2RAD;
    double lam = rec_lon_deg * DEG2RAD;
    double sp = sin(phi), cp = cos(phi);
    double sl = sin(lam), cl = cos(lam);


    // Rotate the Relative vector dx, dy, dz to ENU
    double E = -sl * dx + cl * dy;
    double N = -sp * cl * dx - sp * sl * dy + cp * dz;
    double U = cp * cl * dx + cp * sl * dy + sp * dz;

    // Range, elevation, azimuth
    range_m = sqrt(E * E + N * N + U * U);
    el_deg = asin(U / range_m) * RAD2DEG;
    az_deg = atan2(E, N) * RAD2DEG;
    if (az_deg < 0) az_deg += 360.0;  // make positive
}


// New Today: Pick the best ephemeris for a satellite at time t
// Returns iterators to the ephemeris closest in time to t_gps for this prn
// findBestEphemeris() answers this question:
// "I have hundreds of ephemerides in my vector. Out of all of them, 
// which one belongs to satellite #5 AND has a time-of-ephemeris closest to right now?"

// Why this matters: Remember from Day 3, your BRDC file has 440 ephemerides(about 14 per satellite).
// For each satellite, you have multiple versions recorded throughout the day.You want the one that's
// most accurate for your target time — the one whose toe is closest to t_gps.
// This function searches the entire vector and returns a pointer (address in memory) to the best match.

//const: "Whatever I return, the caller can't modify it"
//GpsEphemeris: The type of thing we're pointing to. 
// *: POINTER — an address (where it is located in memory), not the data itself
// So this function returns a pointer to a GpsEphemeris. NOT a copy.NOT the struct itself. Just its address in memory.
const GpsEphemeris* findBestEphemeris(const vector<GpsEphemeris>& ephs, int prn, double t_gps)
{
    const GpsEphemeris* best = nullptr;      //nullptr means "no pointer yet — empty". We start with nothing found, then update as we discover matches.
    double best_dt = 1e18;                  //will hold the time difference of the best match so far. We want to find the SMALLEST dt, so we start at a deliberately HUGE value.
    for (const auto& e : ephs) {            //read-only reference to each element in the ALL 440 ephemerides. 
        if (e.prn != prn) continue;
        double dt = abs(e.toe - t_gps);
        if (dt > 302400.0) dt = 604800.0 - dt;  // week rollover
        if (dt < best_dt) { best_dt = dt; best = &e; }  //The & operator means "give me the address of". So, '&e' gives the address of e in the memory
    }
    return best;
}

// ── MAIN ─────────────────────────────────────────────────────────────────
int main()
{
    cout << "=== All Visible GPS Satellites ===" << endl << endl;

    // STEP 1: Define receiver location (U of M, Memphis)
    double rec_lat = 35.1186;       // degrees N
    double rec_lon = -89.9387;      // degrees (west = negative)
    double rec_h = 76.0;          // meters above ellipsoid
    cout << "Receiver location: U of M, Memphis" << endl;
    cout << "  Latitude:  " << rec_lat << " deg" << endl;
    cout << "  Longitude: " << rec_lon << " deg" << endl;
    cout << "  Height:    " << rec_h << " m" << endl << endl;


    SatPosition rec = llaToEcef(rec_lat, rec_lon, rec_h);
    cout << "Receiver ECEF: (" << rec.X << ", " << rec.Y << ", " << rec.Z << ")" << endl << endl;


    // STEP 2: Parse navigation file
    string filename = "../../data/BRDC00IGS_R_20240070000_01D_MN.rnx";
    vector<GpsEphemeris> ephs = parseRinexNav(filename);
    if (ephs.empty()) { cerr << "No ephemerides!" << endl; return 1; }
    cout << "Parsed " << ephs.size() << " ephemerides." << endl << endl;

    // STEP 3: Pick a target time — say, 12:00:00 UTC = 43200 s of GPS week
    // (Sunday Jan 7 2024 was the start of GPS week — so 12:00 UTC = 43200 s)
    double t_gps = 43200.0;   // Sunday noon
    cout << "Target time: " << t_gps << " s of GPS week (12:00 UTC)" << endl << endl;

    // STEP 4: For each unique PRN, find best eph and compute position + sky angles
    vector<VisibleSat> visible;
    for (int prn = 1; prn <= 32; prn++)
    {
        const GpsEphemeris* eph = findBestEphemeris(ephs, prn, t_gps);
        if (!eph) continue;  // no ephemeris for this PRN

        SatPosition sp = computeSatPosECEF(*eph, t_gps);  // *eph: Go to the address stored in eph and access the actual GpsEphemeris object.

        VisibleSat vs;
        vs.prn = prn;
        vs.X = sp.X; vs.Y = sp.Y; vs.Z = sp.Z;
        computeElevAzim(sp.X, sp.Y, sp.Z, rec.X, rec.Y, rec.Z,
            rec_lat, rec_lon,
            vs.elevation_deg, vs.azimuth_deg, vs.range_m);

        // Only keep satellites above 5° elevation (typical mask)
        if (vs.elevation_deg > 5.0)
        {
            visible.push_back(vs);
        }

    }

    // STEP 5: Sort by elevation (highest first)
    // Insted of writing:
        // bool compareElevation(const VisibleSat& a, const VisibleSat& b)
        //{
        //    return a.elevation_deg > b.elevation_deg;
        //}
    // and then:
        // sort(visible.begin(), visible.end(), compareElevation);

    // you can write the comparison function directly where it's needed using lambda function:
    sort(visible.begin(), visible.end(), [](const VisibleSat& a, const VisibleSat& b) 
        {
            return a.elevation_deg > b.elevation_deg;
        });

    //In C++, every lambda expression must start with [], even if it doesn't capture anything.
    // lambda structure: 
        // [captures](parameters) 
        // {
        //    body;
        // }



    // STEP 6: Print the sky-view table
    cout << "VISIBLE GPS SATELLITES (elevation > 5 deg):" << endl;
    cout << "==================================================================" << endl;
    cout << setw(4) << "PRN"
        << setw(12) << "Elev (deg)"
        << setw(12) << "Azim (deg)"
        << setw(15) << "Range (km)" << endl;
    cout << "------------------------------------------------------------------" << endl;
    cout << fixed << setprecision(2);
    for (const auto& v : visible) {
        cout << "G" << setw(2) << setfill('0') << v.prn << setfill(' ')
            << setw(13) << v.elevation_deg
            << setw(12) << v.azimuth_deg
            << setw(15) << v.range_m / 1000.0
            << endl;
    }
    cout << "==================================================================" << endl;
    cout << "Total visible satellites: " << visible.size() << endl;

    return 0;

}