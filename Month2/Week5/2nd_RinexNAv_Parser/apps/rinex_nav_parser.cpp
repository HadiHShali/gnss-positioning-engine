// ===========================================================================//
// rinex_nav_parser.cpp - Parse GPS broadcast ephemeris from Rinex 3 Nav file //
// ===========================================================================//

// Reads BRDC file, extract all GPS satellite ephemeris into a vector

// SECTION0: Headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Physical Constants
const double MU_EARTH = 3.986005e14;		// GM constant
const double OMEGA_EARTH = 7.2921151467e-5; // Earth's rotation rate (rad/sec)

// Data structure
struct GpsEphemeris
{
	int prn;		 // Satellite PRN (1-32)
	int year, month, day, hour, minute;
	double second;	 // epoch of clock parameters
	double clk_bias; // SV clock bias (s)
	double clk_drift; // SV clock drift (s/s)
	double clk_drift_rate; // SV clock drift rate(s/s^2)
// Keplerian orbit parameters and corrections
	double IODE;         // Issue of Data, Ephemeris (ephemeris version number)
	double Crs;          // Sine harmonic correction to orbit radius (m)
	double delta_n;      // Mean motion difference from computed value (rad/s)
	double M0;           // Mean anomaly at reference time toe (rad)

	double Cuc;          // Cosine harmonic correction to argument of latitude (rad)
	double e;            // Orbital eccentricity
	double Cus;          // Sine harmonic correction to argument of latitude (rad)
	double sqrt_a;       // Square root of semi-major axis (sqrt(m))

	double toe;          // Reference epoch of ephemeris (s of GPS week)
	double Cic;          // Cosine harmonic correction to inclination angle (rad)
	double Omega0;       // Longitude of ascending node at weekly epoch (rad)
	double Cis;          // Sine harmonic correction to inclination angle (rad)

	double i0;           // Inclination angle at reference time toe (rad)
	double Crc;          // Cosine harmonic correction to orbit radius (m)
	double omega;        // Argument of perigee (rad)
	double Omega_dot;    // Rate of right ascension (rad/s)

	double i_dot;        // Rate of inclination angle (rad/s)
	double L2_codes;     // L2 channel code indicator
	double gps_week;     // GPS week number
	double L2P_flag;     // L2 P-code data flag

	double sv_accuracy;  // User Range Accuracy (URA) index
	double sv_health;    // Satellite health status (0 = healthy)
	double TGD;          // Group Delay Differential (L1-L2 bias) (s)
	double IODC;         // Issue of Data, Clock (clock data version)

	double trans_time;   // Navigation message transmission time (s of GPS week)
	double fit_interval; // Ephemeris fit interval (hours)

};

struct SatPosition
{
	double X, Y, Z;
};

// Helper Function
// Convert Fortran 'D' notation to C++ 'E' notation
string fortranToCpp(const string& s)
{
	string out = s;
	for (char& c : out)
	{
		if (c == 'D' || c == 'd') c = 'E';
	}
	return out;
}

// Helper: safely parse a fixed-width field as double
double parseField(const string& line, size_t start, size_t len = 19)
{
	if (start >= line.length()) return 0.0;
	string field = line.substr(start, len);
	field = fortranToCpp(field);
	try { return stod(field); }
	catch (...) { return 0.0; }
}

// Read 4 numbers from a broadcast orbit line
void readOrbitLine(const string& line, double& v1, double& v2, double& v3, double& v4)
{
	v1 = parseField(line, 4);
	v2 = parseField(line, 23);
	v3 = parseField(line, 42);
	v4 = parseField(line, 61);
}

// Main Parser
vector<GpsEphemeris> parseRinexNav(const string& filename)
{
	vector<GpsEphemeris> ephemerides;
	ifstream file(filename);
	if (!file.is_open())
	{
		cerr << "Error: Can not Open the" << filename << endl;
		return ephemerides;
	}

	string line;
	bool header_done = false;

	// skip the header
	while (getline(file, line))
	{
		if (line.find("END OF HEADER") != string::npos)
		{
			header_done = true;
			break;
		}
	}
	if (!header_done)
	{
		cerr << "Error: No END OF HEADER found" << endl;
		return ephemerides;
	}


	// Step2 : Read records
	while (getline(file, line))
	{
		//skip empty lines
		if (line.length() < 23) continue;

		// GPS records start with 'G'
		if (line[0] != 'G')
		{
			// skip non-GPS satellite for today (R, E, C. etc.)
			// each line is 8 lines - skip 7 more lines
			for (int i = 0; i < 7; i++)
			{
				if (!getline(file, line)) break;//!getline(file, line) becomes true if: End-of-file (EOF) is reached, or A read error occurs.
			}
			continue;
		}

		// Step3: We have a GPS record. Parse all 8 line.
		GpsEphemeris eph = {}; // zero initialize
		// line 0: PRN, Epoch, Clock parameters
		try
		{
			eph.prn = stoi(line.substr(1, 2));
			eph.year = stoi(line.substr(4, 4));
			eph.month = stoi(line.substr(9, 2));
			eph.day = stoi(line.substr(12, 2));
			eph.hour = stoi(line.substr(15, 2));
			eph.minute = stoi(line.substr(18, 2));
			eph.second = parseField(line, 21, 2);
		}
		catch (...) { continue; }

		eph.clk_bias = parseField(line, 23);
		eph.clk_drift = parseField(line, 42);
		eph.clk_drift_rate = parseField(line, 61);


		// Lines 1-7: broadcast orbit parameters (4 numbers each)
		if (!getline(file, line)) break;
		readOrbitLine(line, eph.IODE, eph.Crs, eph.delta_n, eph.M0);

		if (!getline(file, line)) break;
		readOrbitLine(line, eph.Cuc, eph.e, eph.Cus, eph.sqrt_a);

		if (!getline(file, line)) break;
		readOrbitLine(line, eph.toe, eph.Cic, eph.Omega0, eph.Cis);

		if (!getline(file, line)) break;
		readOrbitLine(line, eph.i0, eph.Crc, eph.omega, eph.Omega_dot);

		if (!getline(file, line)) break;
		readOrbitLine(line, eph.i_dot, eph.L2_codes, eph.gps_week, eph.L2P_flag);

		if (!getline(file, line)) break;
		readOrbitLine(line, eph.sv_accuracy, eph.sv_health, eph.TGD, eph.IODC);

		if (!getline(file, line)) break;
		readOrbitLine(line, eph.trans_time, eph.fit_interval,
			*(new double), *(new double));  // last 2 are spare


		// Step 4: Store the complete ephemeris
		ephemerides.push_back(eph); 
	}

	file.close();
	return ephemerides; 
}

// Kepler solver
// ── KEPLER SOLVER (from Day 2) ───────────────────────────────────────────
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

// MAin Function
int main()
{
	cout << "===RINEX Nav Parser + GPS Satellite Position===" << endl << endl;

	string filename = "../../data/BRDC00IGS_R_20240070000_01D_MN.rnx";
	vector<GpsEphemeris> ephs = parseRinexNav(filename);

	cout << "Successfully Parsed " << ephs.size() << " GPS Ephemeris records." << endl << endl;

	if (ephs.empty())
	{
		cerr << "No ephemeris found. Exiting." << endl;
		return 1;
	}

	// Show summary of first 3
	cout << "First 3 ephemerides:" << endl;
	cout << "PRN  toe(s)    sqrt_a(m^0.5)  e        i0(rad)" << endl;
	cout << "-------------------------------------------------------" << endl;
	for (size_t i = 0; i < min((size_t)3, ephs.size()); i++) {
		const auto& e = ephs[i];  //compiler deduces the type (probably Ephemeris).
		cout << "G" << setw(2) << setfill('0') << e.prn << "  "  //setw(2): Set the width of the next output field to 2 characters. setfill('0'): Fill unused spaces with zeros.
			<< fixed << setprecision(1) << setw(7) << e.toe << "  " //fixed: Use fixed-point notation for floating-point numbers. setprecision(1): Show 1 digit after the decimal point.
			<< setprecision(6) << setw(13) << e.sqrt_a << "  "
			<< scientific << setprecision(3) << e.e << "  " //scientific: Switch to scientific notation.
			<< fixed << setprecision(4) << e.i0 << endl;
	}
	cout << endl;

	// Compute position of the first satellite at its own toe
	const GpsEphemeris& eph = ephs[0];
	SatPosition pos = computeSatPosECEF(eph, eph.toe);

	double r = sqrt(pos.X * pos.X + pos.Y * pos.Y + pos.Z * pos.Z);
	cout << "Sample computation PRN " << eph.prn
		<< " at t = " << eph.toe << " s of week:" << endl;
	cout << fixed << setprecision(2);
	cout << "  X = " << pos.X << " m" << endl;
	cout << "  Y = " << pos.Y << " m" << endl;
	cout << "  Z = " << pos.Z << " m" << endl;
	cout << "  |r| = " << r / 1000.0 << " km (expected ~26,560)" << endl;

	return 0;

}