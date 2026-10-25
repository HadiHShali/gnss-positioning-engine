// clk_reader_crctd.cpp - corrected copy of clk_reader.cpp (Month 3 / Week 14 review, 2026-10-24)
// The original stays untouched. Changes are marked [CRCTD].
//
// WHAT WAS WRONG: the reader took EVERY line starting with "AS" as a satellite clock record,
// including HEADER lines. The header lists the reference stations, one per line, starting
// with the 4-character station name:
//     ASCG 30602M004            6121151559 -1563978988  -872615266SOLN STA NAME / NUM
//     ASPA 50503S006           -6100260229  -996502364 -1567977079SOLN STA NAME / NUM
// Ascension Island (ASCG) and Pago Pago (ASPA) both start with "AS". For those lines,
// `ss >> yr` hit 6121151559 (too big for an int), the stream failed, and hr, mi, sec and
// clk_bias were never set (uninitialized = undefined behavior). The result was an extra,
// random-time epoch: 289 CLK epochs instead of 288.
// EFFECT ON RESULTS: none here - the "PRN" of those records was the station number
// (e.g. "30602M004"), which never matches "G01".."G32". But it was luck, not design.
//
// CORRECTIONS:
//   1. Records are read only AFTER "END OF HEADER".
//   2. A record must start with "AS " (AS + space), not just "AS".
//   3. If any field fails to parse, the line is skipped (no uninitialized values).
//   4. A file that cannot be opened is reported (was: silently 0 epochs).
// NOT CHANGED: t_gps is still seconds of DAY (correct for Sunday data only, like the
//   RINEX and SP3 readers).
#include <clk_reader.h>
#include <fstream>
#include <sstream>
#include <iostream>   // [CRCTD] for the error message
using namespace std;

vector<ClkEpoch> parseClk(const string& filename)
{
	vector<ClkEpoch> epochs;

	ifstream file(filename);
	if (!file.is_open()) {                                       // [CRCTD] 4
		cerr << "Cannot open CLK file: " << filename << endl;
		return epochs;
	}

	string line;

	// [CRCTD] 1. skip the header, so station lines like "ASCG ... SOLN STA NAME / NUM" are never read as data
	while (getline(file, line)) {
		if (line.find("END OF HEADER") != string::npos) break;
	}

	map<double, ClkEpoch> by_time; // group AS lines sharing the same epoch
	// groups CLK file records by their timestamp — so all satellites' clock corrections
	// for the same epoch end up bundled together in one ClkEpoch, even though the raw
	// file lists them as separate AS lines, one per satellite.
	// Key type: double — the epoch's GPS time (t_gps)
	// Value type: ClkEpoch — your struct holding t_gps plus the sat_clock map (PRN → clock bias)

    while (getline(file, line)) {
        if (line.compare(0, 3, "AS ") != 0) continue;   // [CRCTD] 2. was substr(0, 2) != "AS"; skips AR (receiver) lines

        istringstream ss(line);
        string tag, prn;
        int yr, mo, dy, hr, mi, nvals;
        double sec, clk_bias;

        //This line reads 10 space-separated pieces of text off the line, one after another,
        // automatically converting each into the correct variable type (string, int, or double) as it goes
        // [CRCTD] 3. the whole read is checked: if any field fails, skip the line
        if (!(ss >> tag >> prn >> yr >> mo >> dy >> hr >> mi >> sec >> nvals >> clk_bias)) continue;
        //Exa: ↑      ↑      ↑    ↑     ↑
        //     AS    G01    2024  1     7     0     0  0.000000   2   -1.234567890123E-04   1.2E-11 (we don't read the last one)

        // clk_bias is ALREADY in seconds -- no unit conversion needed here

        double t_gps = hr * 3600.0 + mi * 60.0 + sec;   // seconds of DAY (Sunday data only)
        by_time[t_gps].t_gps = t_gps;
        by_time[t_gps].sat_clock[prn] = clk_bias;
    }
    for (auto& [t, ep] : by_time) epochs.push_back(ep);
    return epochs;
}
