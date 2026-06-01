//=============================================//
//				nmea_parser.cpp ===============//
//=============================================//

// SECTION0 : Headers
#include <iostream>
#include <fstream>  // for reading input files : ifstrem
#include <sstream>  // for string spliting : stringstream
#include <string>   // for std::string
#include <cmath>    // for floor
#include <vector>   // for vector to hold fields 

using namespace std;

// Helper Function 1
// Convert NMEA lat/lon to decimal degrees
// NMEA format is DDMM.mmm (degree  + minutes). We need decimal degrees
double nmea_to_decimal(double nmea_coord, char direction) {
	double degrees = floor(nmea_coord / 100.0);
	double minutes = nmea_coord - degrees * 100.0;
	double decimal = degrees + minutes / 60.0;

	// South and West are negative
	if (direction == 'S' || direction == 'W') {
		decimal = -decimal;
	}
	return decimal;
}


// Helper Function 2
// Split a string by delimeter into a vector
// vector <string> : The function RETURNS a vector (list) of strings
vector <string> split(const string& line, char delim) // string& : don't make a copy of the input, just point to it. const: and I won't modify the original
{
	vector <string> fields; //We create an EMPTY vector named 'fields' to hold the pieces we'll cut out. A dynamic list: it can grow as we move on
	stringstream ss(line); // Wrap the string in a stream. It makes the line string behave like a stream of characters that you can read from piece by piece
	
	// declare a temporary variable with a type of string (a "bucket" we'll reuse.)
	string field;

	//  getline(ss, field, delim) does this:
	//		Reads characters from ss(the stream)
	//		Stops when it hits delim(the comma)
	//		Stores what it read into field
	//		Moves the cursor past the comma so next call starts after it
	//		Returns true if it read something, false if the stream is empty
	while (getline(ss, field, delim))  // Keeps repeating as long as 'getline' returns true (i.e., as long as there's more to read)
	{
		fields.push_back(field);
	}
	return fields;
}


// SECTION 1: Main program
int main()
{
	// open the NMEA file
	ifstream nmea_file("../../data/sample_nmea.txt");
	if (!nmea_file.is_open())
	{
		cerr << "ERROR: Could not open sample_nmea.txt" << endl;
		return 1;  //return 1 (or any non-zero number) tells the operating system: "The program failed."
				   //return 0 tells the operating system: "The program succeeded."
	}
	cout << "=== NMEA Parser ===" << endl;
	cout << "Time      Latitude     Longitude    Alt(m)  Sats  HDOP" << endl;
	cout << "-------------------------------------------------------------" << endl;

	// Read the file line by line
	string line;
	int gpgga_count = 0;

	while (getline(nmea_file, line))
	{
		// We only care about the $GPGGA sentences for today
		if (line.substr(0, 6) != "$GPGGA")
		{
			continue; // skip other sentences type
		}

		// split the line by comma
		vector<string> fields = split(line, ',');

		// GPGGA has 15 fields 
		if (fields.size() < 15)
		{
			cerr << "WARNING: malformed GPGGA sentene, skipping." << endl;
			continue;
		}

		// Extract the fields we care about
				// Extract the fields we care about
		string time_utc = fields[1];           // "123519"
		double lat_nmea = stod(fields[2]);     // 4807.038 ,stod: string to double
		char   lat_dir = fields[3][0];        // 'N'
		double lon_nmea = stod(fields[4]);     // 01131.000
		char   lon_dir = fields[5][0];        // 'E'
		int    fix_qual = stoi(fields[6]);     // 1
		int    n_sats = stoi(fields[7]);     // 8
		double hdop = stod(fields[8]);     // 0.9
		double altitude = stod(fields[9]);     // 545.4

		// Convert to decimal degrees
		double lat_deg = nmea_to_decimal(lat_nmea, lat_dir);
		double lon_deg = nmea_to_decimal(lon_nmea, lon_dir);

		// Print the parsed position
		cout << time_utc
			<< "   " << lat_deg
			<< "   " << lon_deg
			<< "   " << altitude
			<< "   " << n_sats
			<< "   " << hdop << endl;

		gpgga_count++;
	}
	nmea_file.close();
	cout << endl << "Parsed " << gpgga_count << " GPGGA sentences." << endl;
	return 0;
}