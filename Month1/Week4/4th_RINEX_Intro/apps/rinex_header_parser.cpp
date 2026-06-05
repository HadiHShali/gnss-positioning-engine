// ==========================================================================//
//            rinex_header_parser.cpp - Read a RINEX 3 OBS file header       //
// ==========================================================================//

// SECTION0: Headers-----------------------------------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <windows.h>
using namespace std;

// SECTION1: Data structure - holds all header info-----------------------------
struct RinexHeader
{
	string version;
	string file_type;
	string station_name;
	string observer;
	string agency;
	string receiver_type;
	string antenna_type;
	double approx_x = 0.0;
	double approx_y = 0.0;
	double approx_z = 0.0;
	int n_obs_types = 0;
	vector<string> obs_types;
	double interval_sec = 0.0;
	string time_first_obs;
};

// SECTION2: Helper Functions------------------------------------------------------
// Helper 1: trim
// This function removes leading (front) and trailing (back) whitespace from a string, but keeps internal spaces intact.
// trim("   he llo   ")        →  "he llo"
string trim(const string& s)
{
	// s.find_first_not_of(" \t"): This method searches the string and returns the position of the first character that is NOT in the given set.
	size_t start = s.find_first_not_of(" \t");
	//String:  "   CERI   "
	//Position :0123456789
	//         ''''CERI'''''''
	//             ↑
	//             Position 3 — first non - whitespace character!
	// So find_first_not_of(" \t"): Find first character that is NOT one of the "space and tab".  returns 3

	//size_t: is an unsigned integer type used for sizes/positions in C++. Think of it as a special int. 

	if (start == string::npos) return "";  //string::npos  >> It's a special value meaning "not found."


	size_t end = s.find_last_not_of(" \t"); //Find the last non-whitespace character. Same idea as above, but searches from the END backwards.
	//String:  "   CERI   "
	//Position :0123456789
	//         ''''CERI'''''''
	//                ↑
	//                Position 6 — first non - whitespace character!
	// So find_last_not_of(" \t"): Find last character that is NOT one of the "space and tab".  returns 6
	
	
	return s.substr(start, end - start + 1); 
}

// Helper2 : safely extract a substring
string safe_substr(const string& s, size_t start, size_t len)
{
	if (start >= s.length()) return "";
	return s.substr(start, len);
}

// SECTION3: Main Parser Function ------------------------------------------------------
bool parse_header(const string& filename, RinexHeader& header)
{
	ifstream file(filename);
	if (!file.is_open())
	{
		cerr << "Error: Can Not Open the " << filename << endl;
		return false;
	}

	string line;
	while (getline(file, line))
	{
		// Header lines have data in column 1-60 and label in column 61-80
		// Extract the label (in practice it is going to be from column 60-79 because the indexing starts from zero)
		string label = safe_substr(line, 60, 20);
		label = trim(label);
		string data = safe_substr(line, 0, 60);

		// Check each known header type
		if (label == "RINEX VERSION / TYPE")
		{
			header.version = trim(safe_substr(line, 0, 9));
			header.file_type = trim(safe_substr(line, 20, 20));
		}

		else if (label == "MARKER NAME")
		{
			header.station_name = trim(data);
		}
		else if (label == "OBSERVER / AGENCY") {
			header.observer = trim(safe_substr(line, 0, 20));
			header.agency = trim(safe_substr(line, 20, 40));
		}
		else if (label == "REC # / TYPE / VERS") {
			header.receiver_type = trim(safe_substr(line, 20, 20));
		}
		else if (label == "ANT # / TYPE") {
			header.antenna_type = trim(safe_substr(line, 20, 20));
		}
		else if (label == "APPROX POSITION XYZ")
		{
			//  -252840.2412 -4870321.3215  4106503.1234                  
			// three 14.4 floats: cols 00-13, 14-27, 28-41
			try
			{
				header.approx_x = stod(safe_substr(line, 0, 14));
				header.approx_y = stod(safe_substr(line, 14, 14));
				header.approx_z = stod(safe_substr(line, 28, 14));
			}
			catch (...)  // The compiler ignores everything between /* and */, no matter how many lines
			{
				/* keep zeros */
			}
		}
		else if (label == "SYS / # / OBS TYPES") {
			// First char = system, then count, then obs codes
			try {
				header.n_obs_types = stoi(safe_substr(line, 3, 3));
			}
			catch (...) {}
			// Parse obs codes (3-char each, starting at col 7)
			stringstream ss(safe_substr(line, 7, 60));
			string code;
			while (ss >> code) {  //ss >> code reads one word at a time from a stringstream, automatically splitting on whitespace (spaces, tabs, newlines).
				header.obs_types.push_back(code);
			}
		}
		else if (label == "INTERVAL") {
			try {
				header.interval_sec = stod(trim(data));
			}
			catch (...) {}
		}
		else if (label == "TIME OF FIRST OBS") {
			header.time_first_obs = trim(data);
		}
		else if (label == "END OF HEADER") {
			file.close();
			return true;  // success!
		}
	}
	cerr << "WARNING: Reached end of file without finding END OF HEADER" << endl;
	file.close();
	return false;

}

// SECTION4: PRINT SUMMARY — display the parsed header nicely 
void print_header(const RinexHeader& h) {
	cout << "=======================================================" << endl;
	cout << "  RINEX HEADER SUMMARY" << endl;
	cout << "=======================================================" << endl;
	cout << "Version:        " << h.version << " (" << h.file_type << ")" << endl;
	cout << "Station name:   " << h.station_name << endl;
	cout << "Observer:       " << h.observer << endl;
	cout << "Agency:         " << h.agency << endl;
	cout << "Receiver:       " << h.receiver_type << endl;
	cout << "Antenna:        " << h.antenna_type << endl;
	cout << endl;
	cout << "ECEF Position:" << endl;
	cout << "  X = " << h.approx_x << " m" << endl;
	cout << "  Y = " << h.approx_y << " m" << endl;
	cout << "  Z = " << h.approx_z << " m" << endl;
	cout << endl;
	cout << "Observation types (" << h.n_obs_types << " expected, "
		<< h.obs_types.size() << " parsed):" << endl;
	for (const string& code : h.obs_types) {
		cout << "  - " << code << endl;
	}
	cout << endl;
	cout << "Sample interval: " << h.interval_sec << " sec" << endl;
	cout << "First epoch:     " << h.time_first_obs << endl;
	cout << "=======================================================" << endl;
}


// MAin function
int main() {
	cout << "=== RINEX 3 Header Parser ===" << endl << endl;

	RinexHeader header;
	bool ok = parse_header("../../data/sample.24o", header);

	if (ok) {
		print_header(header);
	}
	else {
		cerr << "Failed to parse header." << endl;
		return 1;
	}

	return 0;
}
