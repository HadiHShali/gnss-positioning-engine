// ====================================================================//
// nmea_parser_v2.cpp, handles GPGGA + GPRMC sentences, and export CSV //
// ====================================================================//

// SECTION0: Headers
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

// SECTION 1: HElper Functions
// Helper 1: Convert NMEA DDMM.mmmm to decimal degrees
// ── HELPER: Convert NMEA DDMM.mmmm to decimal degrees ────────────────────
double nmea_to_decimal(double nmea_coord, char direction) {
	double degrees = floor(nmea_coord / 100.0);
	double minutes = nmea_coord - degrees * 100.0;
	double decimal = degrees + minutes / 60.0;
	if (direction == 'S' || direction == 'W') decimal = -decimal;
	return decimal;
}


// Helper 2: Split a string by delimiter into a vector
vector<string> split(const string& line, char delim)
{
	vector<string> fields; // output list of strings: empty vector that needs to be filled: []
	string field_dummy;
	stringstream ss(line);
	// reads charachter from ss, stops when it hit the delim, put the value in field_dummy variable, and then return true if it reads something
	while (getline(ss, field_dummy, delim)) 
	{
		fields.push_back(field_dummy);
	}
	return fields;
}

// Helper 3: safe string to double conversion - 
double safe_stod(const string& s)
{
//return 0 if string is empty >> that's the fast path. 
// We handle those exceptions that are common and let the program handles others. It makes the program fatster
	if (s.empty()) return 0.0;
	try { return stod(s); }
	catch (...) { return 0.0; }//The ... means "catch ANY type of exception"
}

// Section 3: Main Program
int main()
{
	// open the input file
	ifstream input("../../data/walk.nmea");
	if (!input.is_open())
	{
		cerr << "Error: The input file can not be open";
		return 1;
	}

	// open the output csv file
	ofstream output("../../data/gps_track.csv");
	if (!output.is_open())
	{
		cerr << "Error: Can not creat the output file";
		return 1;
	}

	// write CSV header
	output << "time,date,lat,lon,alt,sats,hdop,speed_ms,course" << endl;

	cout << "=== NMEA Parser v2 — Real GPS Data ==="  << endl << endl;

	// Counters for summary
	int total_lines = 0;
	int gga_count = 0;
	int rmc_count = 0;
	int skipped = 0;

	// Now we have to combine GGA and RMC sentences at the same time
	// Variables to combine
	// GGA gives alt, sats, hdop , RMC gives date, speed, course
	string last_time = "";
	string last_date = "";
	double last_lat = 0.0, last_lon = 0.0, last_alt = 0.0;
	int last_sats = 0;
	double last_hdop = 0.0, last_speed = 0.0, last_course = 0.0;
	bool have_gga = false, have_rmc = false;

	// loop
	string line;
	while (getline(input, line))
	{
		total_lines++;

		// Skip empty or too short line
		if (line.length() < 7) { skipped++; continue; }

		// Must start with $ sign
		if (line[0] != '$') { skipped++; continue; }

		// identify sentence type
		string id = line.substr(0, 6);
		vector<string> f = split(line, ',');

		// GPGGA: position + altitude + satellites
		if (id == "$GPGGA" && f.size()>= 10)
		{
			//Checking the fix quaklity	
			if (safe_stod(f[6]) == 0) { skipped++; continue; }

			last_time = f[1]; //UTC time
			last_lat  = nmea_to_decimal(safe_stod(f[2]), f[3].empty() ? 'N' : f[3][0]); //IF f[3] is empty, use 'N'. OTHERWISE, use the first character of f[3].
			last_lon  = nmea_to_decimal(safe_stod(f[4]), f[5].empty() ? 'E' : f[5][0]);
			last_sats = int(safe_stod(f[7]));
			last_hdop = safe_stod(f[8]);
			last_alt  = safe_stod(f[9]);
			have_gga  = true;
			gga_count++;
		}
		else if (id == "$GPRMC" && f.size() >= 10)
		{
			// Check the status (field 2): A=active, V=void
			if (f[2].empty() || f[2][0] != 'A') { skipped++; continue; }
			
			last_time = f[1];
			last_lat = nmea_to_decimal(safe_stod(f[3]), f[4].empty() ? 'N' : f[4][0]);
			last_lon = nmea_to_decimal(safe_stod(f[5]), f[6].empty() ? 'E' : f[6][0]);
			// Speed is in knots — convert to m/s (1 knot = 0.5144 m/s)
			last_speed = safe_stod(f[7]) * 0.5144;
			last_course = safe_stod(f[8]);
			last_date = f[9];
			have_rmc = true;
			rmc_count++;
		}
		// Skip everything else (GPGSV, GPGSA, GPVTG, etc.)
		else 
		{
			skipped++;
			continue;
		}

		// Write a CSV row whenever we have both GGA and RMC for this epoch
		if (have_gga && have_rmc)
		{
			output << last_time << ","
				<< last_date << ","
				<< last_lat << ","
				<< last_lon << ","
				<< last_alt << ","
				<< last_sats << ","
				<< last_hdop << ","
				<< last_speed << ","
				<< last_course << endl;
		// Reset for the next epoch
			have_gga = false;
			have_rmc = false;
		}
	}

	input.close();
	output.close();

	// Print summary
	cout << "Total lines read:       " << total_lines << endl;
	cout << "GPGGA sentences parsed: " << gga_count << endl;
	cout << "GPRMC sentences parsed: " << rmc_count << endl;
	cout << "Skipped lines:          " << skipped << endl;
	cout << endl << "Output written to: data/gps_track.csv" << endl;
	return 0;

}