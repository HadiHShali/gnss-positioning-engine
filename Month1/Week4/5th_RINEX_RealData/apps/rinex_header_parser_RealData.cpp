// ==========================================================================//
//            rinex_header_parser.cpp - Read a RINEX 3 OBS file header       //
// ==========================================================================//

// SECTION0: Headers-----------------------------------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
using namespace std;

// SECTION1: Data structure - holds all header info-----------------------------
struct RinexHeader {
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
	double interval_sec = 0.0;
	string time_first_obs;
	string time_last_obs;
	map<char, vector<string>> obs_types;  // per-system
};


// SECTION2: Helper Functions------------------------------------------------------
// Helper 1: trim
// This function removes leading (front) and trailing (back) whitespace from a string, but keeps internal spaces intact.
/// ── HELPERS ──────────────────────────────────────────────────────────────
string trim(const string& s) {
    size_t start = s.find_first_not_of(" \t");
    if (start == string::npos) return "";
    size_t end = s.find_last_not_of(" \t");
    return s.substr(start, end - start + 1);
}

string safe_substr(const string& s, size_t start, size_t len) {
    if (start >= s.length()) return "";
    return s.substr(start, len);
}

// Map constellation letter to friendly name
string constellation_name(char c) {
	switch (c) {
	case 'G': return "GPS (USA)";
	case 'R': return "GLONASS (Russia)";
	case 'E': return "Galileo (Europe)";
	case 'C': return "BeiDou (China)";
	case 'J': return "QZSS (Japan)";
	case 'S': return "SBAS";
	case 'I': return "IRNSS/NavIC (India)";
	default:  return "Unknown";
	}
}


// SECTION3: Main Parser Function ------------------------------------------------------
// SECTION3: Main Parser Function -------------------------------------------
bool parse_header(const string& filename, RinexHeader& header)
{
    ifstream file(filename);

    if (!file.is_open())
    {
        cerr << "ERROR: Cannot open file: " << filename << endl;
        return false;
    }

    char current_sys = ' ';
    string line;

    while (getline(file, line))
    {
        string label = trim(safe_substr(line, 60, 20));
        string data = safe_substr(line, 0, 60);

        if (label == "RINEX VERSION / TYPE")
        {
            header.version = trim(safe_substr(line, 0, 9));
            header.file_type = trim(safe_substr(line, 20, 20));
        }

        else if (label == "MARKER NAME")
        {
            header.station_name = trim(data);
        }

        else if (label == "OBSERVER / AGENCY")
        {
            header.observer = trim(safe_substr(line, 0, 20));
            header.agency = trim(safe_substr(line, 20, 40));
        }

        else if (label == "REC # / TYPE / VERS")
        {
            header.receiver_type = trim(safe_substr(line, 20, 20));
        }

        else if (label == "ANT # / TYPE")
        {
            header.antenna_type = trim(safe_substr(line, 20, 20));
        }

        else if (label == "APPROX POSITION XYZ")
        {
            try
            {
                header.approx_x = stod(safe_substr(line, 0, 14));
                header.approx_y = stod(safe_substr(line, 14, 14));
                header.approx_z = stod(safe_substr(line, 28, 14));
            }
            catch (...)
            {
            }
        }

        else if (label == "SYS / # / OBS TYPES")
        {
            // New system line or continuation line
            char first_char = line.empty() ? ' ' : line[0];

            if (first_char != ' ')
            {
                current_sys = first_char;
            }

            // Observation codes begin around column 7
            stringstream ss(safe_substr(line, 6, 54));

            string token;
            while (ss >> token)
            {
                // Observation codes are 3 characters:
                // C1C, L1C, D1C, S1C, etc.
                if (token.length() == 3)
                {
                    header.obs_types[current_sys].push_back(token);
                }
            }
        }

        else if (label == "INTERVAL")
        {
            try
            {
                header.interval_sec = stod(trim(data));
            }
            catch (...)
            {
            }
        }

        else if (label == "TIME OF FIRST OBS")
        {
            header.time_first_obs = trim(data);
        }

        else if (label == "TIME OF LAST OBS")
        {
            header.time_last_obs = trim(data);
        }

        else if (label == "END OF HEADER")
        {
            file.close();
            return true;
        }
    }

    cerr << "WARNING: Reached EOF before END OF HEADER" << endl;
    file.close();
    return false;
}
// ── PRINT SUMMARY ────────────────────────────────────────────────────────
void print_header(const RinexHeader& h) {
    cout << "=======================================================" << endl;
    cout << "  RINEX HEADER SUMMARY (Multi-GNSS)" << endl;
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
    cout << "Sample interval: " << h.interval_sec << " sec" << endl;
    cout << "First epoch:     " << h.time_first_obs << endl;
    cout << "Last epoch:      " << h.time_last_obs << endl;
    cout << endl;
    cout << "------------- Constellations Observed -----------" << endl;
    for (const auto& [sys, codes] : h.obs_types) {
        cout << "  [" << sys << "] " << constellation_name(sys)
            << " - " << codes.size() << " observables: ";
        for (size_t i = 0; i < codes.size(); ++i) {
            cout << codes[i];
            if (i + 1 < codes.size()) cout << ", ";
        }
        cout << endl;
    }
    cout << "----------------------------------------------------" << endl;
}

// ── MAIN ─────────────────────────────────────────────────────────────────
int main() {
    cout << "=== RINEX 3 Multi-GNSS Header Parser ===" << endl << endl;

    RinexHeader header;
    string filename = "../../data/BILL00USA_R_20240070000_01D_15S_MO.crx/BILL00USA_R_20240070000_01D_15S_MO.rnx";

    bool ok = parse_header(filename, header);
    if (ok) {
        print_header(header);
    }
    else {
        cerr << "Failed to parse header." << endl;
        return 1;
    }
    return 0;
}
