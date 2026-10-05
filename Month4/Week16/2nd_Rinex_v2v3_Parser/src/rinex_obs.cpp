// rinex_obs.cpp: version aware RINEX observation reader (2.xx, 3.xx)
//
// Design: the the version only changes three things:
//          a) how the header lists observation types
//			b) how the epoch line looks
//			c) how one satellite's observations are laid out on the page

// Everything else is shared: each satellite's data is turned into one string of 16-char slots (14 values + 1 LLI + 1 signal strength),
// and the same readSlot() pulls values out of it for both versions.

#include "rinex_obs.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

//===========================================================
// small helpers
// =========================================================
// getline that also strips a trailing '\r' (files with Windows line ending)
// static: make it private to this .cpp file
static bool getLine(istream& in, string& line) {
	if (!getline(in, line)) return false;
	if (!line.empty() && line.back() == '\r') line.pop_back();
	return true;
}

static string trim(const string& s) {
	size_t a = s.find_first_not_of(" \t");      // index of the first non-blank character
	if (a == string::npos) return "";           // nothing but blanks -> empty string
	size_t b = s.find_last_not_of(" \t");       // index of the last non-blank character
	return s.substr(a, b - a + 1);              // the part from a to b, inclusive
}

// Fixed-width field. RINEX allows trailing blanks to be stripped, so a short
// line is padded with blanks instead of throwing.
static string field(const string& line, size_t pos, size_t len) {
	if (pos >= line.size()) return string(len, ' ');   // case 1: field starts past the end
	string f = line.substr(pos, len);                   // case 2: normal cut
	if (f.size() < len) f.append(len - f.size(), ' '); // case 3: line ends inside the field
	return f;
}

static int toInt(const string& s, int fallback = 0) {
	string t = trim(s);
	if (t.empty()) return fallback;
	try { return stoi(t); }
	catch (...) { return fallback; }
}

static double toDouble(const string& s, double fallback = 0.0) {
	string t = trim(s);
	if (t.empty()) return fallback;
	try { return stod(t); }
	catch (...) { return fallback; }
}

// Header label sits in columns 61-80
static string headerLabel(const string& line) {
	return (line.size() > 60) ? trim(line.substr(60)) : "";
}

// ============================================================
// GPS time (needs to review)
// ============================================================

// Number of days from 1970-01-01 to the date y-m-d (Gregorian calendar).
// Trick: treat March as the first month of the year, so February (the only
// month whose length changes) is the LAST month and the leap day falls at the end.
static long daysFromCivil(int y, int m, int d) {
	y -= (m <= 2);                       // Jan/Feb belong to the previous "March-year" (true = 1)
	const long era = (y >= 0 ? y : y - 399) / 400;   // which 400-year cycle (the -399 rounds negative years down)
	const long yoe = y - era * 400;      // year within that cycle: 0..399
	const long doy = (153 * (m + (m > 2 ? -3 : 9)) + 2) / 5 + d - 1;   // day of year, Mar 1 = 0:
	//   m + (m > 2 ? -3 : 9) -> Mar=0, Apr=1, ..., Dec=9, Jan=10, Feb=11
	//   (153*mp + 2) / 5     -> days before that month (0, 31, 61, ..., 337)
	//   + d - 1              -> day of the month, counted from 0
	const long doe = yoe * 365 + yoe / 4 - yoe / 100 + doy;   // day within the cycle: 0..146096
	//   365 per year, +1 leap day every 4 years, -1 every 100
	//   (the +1 every 400 is never needed, since yoe < 400)
	return era * 146097 + doe - 719468;  // 146097 days per 400 years; -719468 makes 1970-01-01 = day 0
}

// Fill ep.gps_week and ep.t_gps from the calendar date and time already stored in ep.
static void toGpsTime(ObsEpoch& ep) {
	const long days = daysFromCivil(ep.year, ep.month, ep.day)   // days from 1970-01-01 to this date
		- daysFromCivil(1980, 1, 6);                 // minus days to the GPS epoch = days since GPS start
	ep.gps_week = static_cast<int>(days / 7);                    // whole weeks since 1980-01-06 (explicit long -> int)
	ep.t_gps = (days % 7) * 86400.0                              // day of the week (Sunday = 0) in seconds
		+ ep.hour * 3600.0                                  // + hours in seconds
		+ ep.minute * 60.0                                  // + minutes in seconds
		+ ep.second;                                        // + seconds (can be fractional)
}


// ============================================================
// Version detection (line 1: F9.2 version, file type in column 21)
// ============================================================

// Hatanaka-compressed files carry "CRINEX VERS   / TYPE" on line 1 instead
static bool isHatanaka(const string& line) {
	return headerLabel(line).rfind("CRINEX", 0) == 0;       // true if the label starts with "CRINEX"
}

double detectRinexVersion(const string& filename, char* file_type) {   // public: declared in rinex_obs.h
	ifstream f(filename);                                    // open the file
	string line;                                             // will hold line 1
	if (!f.is_open()) {                                      // wrong path or name
		cerr << "ERROR: cannot open " << filename << endl;
		return 0.0;
	}
	if (!getLine(f, line)) {                                 // file exists but has no lines
		cerr << "ERROR: " << filename << " is empty" << endl;
		return 0.0;
	}
	if (!line.empty() && line[0] == '\x1f') {               // .gz and .Z files both start with byte 0x1F
		cerr << "ERROR: " << filename << " is still compressed (.gz/.Z) -- run: gzip -d <file>" << endl;
		return 0.0;
	}
	if (isHatanaka(line)) {                                  // the case you just met with bill0070.24d
		cerr << "ERROR: " << filename << " is Hatanaka-compressed (CRINEX) -- run: CRX2RNX " << filename << endl;
		return 0.0;
	}
	if (headerLabel(line) != "RINEX VERSION / TYPE") {       // anything else that isn't RINEX
		cerr << "ERROR: " << filename << ": line 1 is not 'RINEX VERSION / TYPE' -- not a RINEX file?" << endl;
		return 0.0;
	}
	if (file_type) *file_type = field(line, 20, 1)[0];       // column 21: 'O', 'N', ...
	return toDouble(field(line, 0, 9));                      // columns 1-9: "     3.04" -> 3.04
}

// ============================================================
// Header
// ============================================================
static bool readHeader(istream& file, RinexObsHeader& h) {        // fills h; returns false on any problem
    string line;                                                  // current header line

    // ---- Line 1: version and file type ----
    if (!getLine(file, line) || headerLabel(line) != "RINEX VERSION / TYPE") {
        cerr << "ERROR: line 1 is not 'RINEX VERSION / TYPE'" << endl;   // safety net (detectRinexVersion checks first)
        return false;
    }
    h.version = toDouble(field(line, 0, 9));                      // "     3.04" -> 3.04
    h.major = static_cast<int>(h.version);                        // 3.04 -> 3   (2.11 -> 2)
    const char type = field(line, 20, 1)[0];                      // column 21: 'O' for observation files
    if (type != 'O') {                                            // e.g. a nav file passed by mistake
        cerr << "ERROR: not an observation file (type '" << type << "')" << endl;
        return false;
    }
    if (h.major < 2 || h.major > 4) {                             // 2.xx and 3.xx (4.xx obs uses the 3.xx layout)
        cerr << "ERROR: RINEX version " << h.version << " not supported" << endl;
        return false;
    }

    char cur_sys = ' ';      // v3: which system the current OBS TYPES block belongs to ('G', 'R', 'E', ...)
    int  n_gps_types = 0;    // how many GPS types the header announces

    // ---- Remaining header lines, until END OF HEADER ----
    while (getLine(file, line)) {
        const string label = headerLabel(line);                   // columns 61-80, trimmed

        if (label == "END OF HEADER") return true;                // success: stream now sits at the first epoch

        if (label == "MARKER NAME") {
            h.marker_name = trim(field(line, 0, 60));             // e.g. "BILL"
        }
        else if (label == "APPROX POSITION XYZ") {                // 3F14.4: three 14-character fields
            h.approx_x = toDouble(field(line, 0, 14));
            h.approx_y = toDouble(field(line, 14, 14));
            h.approx_z = toDouble(field(line, 28, 14));
        }
        else if (label == "INTERVAL") {                           // F10.3, e.g. "    15.000"
            h.interval = toDouble(field(line, 0, 10));
        }
        else if (h.major == 2 && label == "# / TYPES OF OBSERV") {
            // ---------------------------------------------------------------
            // RINEX 2 observation types
            // ---------------------------------------------------------------
            // * ONE list for the whole file. GPS, GLONASS, Galileo and SBAS
            //   satellites all use the same slots; a slot a satellite does not
            //   have is left blank in its data. So for v2, gps_obs_types is
            //   simply "the list" (GPS satellites use it like everyone else).
            // * Codes are 2 characters: L1, L2 (phase), C1, P1, P2 (code),
            //   S1, S2 (signal strength), D1, D2 (Doppler), ...
            // * At most 9 codes fit on one line; more codes continue on the
            //   next line(s) with the same label.
            //
            // Format (RINEX 2.11 spec, Fortran notation):
            //   first line:        I6, 9(4X,A2)
            //   continuation line: 6X, 9(4X,A2)
            //     I6 = integer, right-aligned in 6 columns  (the count)
            //     6X = 6 blank columns                      (no count)
            //     4X = 4 blank columns
            //     A2 = 2 characters                         (one code)
            //   Each slot = 4 blanks + 2-char code = 6 columns.
            //   6 + 9 slots x 6 = 60 columns, then the label from index 60.
            //
            // Example (BILL, bill0070.24o: 20 types -> 9 + 9 + 2 on 3 lines):
            //   index:  0     6   10    16    22    28    34    40    46    52    58
            //          "    20    L1    L2    C1    P2    P1    S1    S2    C2    L5# / TYPES OF OBSERV"
            //          "          C5    S5    L6    C6    S6    L7    C7    S7    L8# / TYPES OF OBSERV"
            //          "          C8    S8                                          # / TYPES OF OBSERV"
            //   (L5 and L8 end at index 59, so they look glued to the label.)
            //
            // Where code i (i = 0..8) starts on ANY of these lines:
            //   6 (count field) + 4 (blanks) + 6*i  =  10 + 6*i
            //   i = 0 -> 10 (L1),  i = 1 -> 16 (L2),  ...,  i = 8 -> 58 (L5)
            //
            // First line vs continuation line:
            //   the count field (index 0-5) is "    20" on the first line and
            //   blank on continuations, so toInt(..., -1) returns 20 or -1.
            //   Only the first line sets the count and starts a fresh list.
            // ---------------------------------------------------------------
            const int n = toInt(field(line, 0, 6), -1);      // 20 on line 1, -1 on continuation lines
            if (n > 0) {                                     // first line of the list:
                n_gps_types = n;                             //   how many codes to collect in total
                h.gps_obs_types.clear();                     //   start a fresh list
            }
            // Read up to 9 slots, but stop once n_gps_types codes are collected,
            // so the empty slots on the last line are never read as codes.
            for (int i = 0; i < 9 && (int)h.gps_obs_types.size() < n_gps_types; i++) {
                const string code = trim(field(line, 10 + 6 * i, 2));   // 2 chars at 10 + 6i
                if (!code.empty()) h.gps_obs_types.push_back(code);     // keep it, in file order
            }
        }
        else if (h.major >= 3 && label == "SYS / # / OBS TYPES") {
            // ---------------------------------------------------------------
            // RINEX 3 observation types
            // ---------------------------------------------------------------
            // * One list PER SYSTEM: a block for G (GPS), then blocks for
            //   R (GLONASS), E (Galileo), C (BeiDou), ... We keep only G.
            // * Codes are 3 characters: <type><band><attribute>
            //     type:      C = code, L = phase, D = Doppler, S = signal strength
            //     band:      1 = L1, 2 = L2, 5 = L5
            //     attribute: tracking mode, e.g. C = C/A, W = semi-codeless P(Y),
            //                L = L2C(L), Q = L5 pilot
            //   e.g. C1W = L1 P(Y) code, L2W = L2 P(Y) phase, L2L = L2C phase
            // * At most 13 codes fit on one line; more codes continue on the
            //   next line(s) with the same label.
            //
            // Format (RINEX 3 spec, Fortran notation):
            //   first line of a block: A1, 2X, I3, 13(1X,A3)
            //   continuation line:     6X, 13(1X,A3)
            //     A1 = 1 character (system letter), 2X = 2 blanks,
            //     I3 = integer in 3 columns (the count), 1X = 1 blank,
            //     A3 = 3 characters (one code)
            //   Each slot = 1 blank + 3-char code = 4 columns.
            //   6 + 13 slots x 4 = 58 columns, then the label from index 60.
            //
            // Example (BILL v3: 22 GPS types -> 13 + 9 on 2 lines, then Galileo):
            //   index:  0  3   7   11  15  19  23  27  31  35  39  43  47  51  55
            //          "G   22 C1C L1C D1C S1C C1W S1W C2W L2W D2W S2W C2L L2L D2L  SYS / # / OBS TYPES"
            //          "       S2L C5Q L5Q D5Q S5Q C1L L1L D1L S1L                  SYS / # / OBS TYPES"
            //          "E   ..."   <- next system's block (skipped)
            //
            // Where code i (i = 0..12) starts on ANY of these lines:
            //   1 (letter) + 2 (blanks) + 3 (count) + 1 (blank) + 4*i  =  7 + 4*i
            //   i = 1 -> 11 (L1C),  i = 4 -> 23 (C1W),  i = 6 -> 31 (C2W),  i = 7 -> 35 (L2W)
            //   -> the Week 15 indices L1C = 1, C1W = 4, C2W = 6, L2W = 7
            //
            // Which system does a line belong to?
            //   a letter in column 1 (index 0) starts a new block -> store it in cur_sys;
            //   a blank in column 1 means "continuation of the current block".
            //   cur_sys survives from one loop pass to the next, which is why it is
            //   declared before the while loop.
            // ---------------------------------------------------------------
            if (line[0] != ' ') {                            // a new system block starts here
                cur_sys = line[0];                           //   remember it for its continuation lines
                if (cur_sys == 'G') {                        //   GPS block:
                    n_gps_types = toInt(field(line, 3, 3));  //     count at index 3-5 ("22")
                    h.gps_obs_types.clear();                 //     start a fresh list
                }
            }
            if (cur_sys != 'G') continue;                    // R, E, C, ... lines: skip to the next header line
            // Read up to 13 slots, stopping once n_gps_types codes are collected.
            for (int i = 0; i < 13 && (int)h.gps_obs_types.size() < n_gps_types; i++) {
                const string code = trim(field(line, 7 + 4 * i, 3));   // 3 chars at 7 + 4i
                if (!code.empty()) h.gps_obs_types.push_back(code);    // keep it, in file order
            }
        }
        // every other header line (COMMENT, ANT # / TYPE, ...) is ignored
    }
    cerr << "ERROR: no END OF HEADER found" << endl;              // file ended inside the header
    return false;
}

// ============================================================
// Observable selection
// ============================================================
//
// The header reader gave us the GPS observation-type list in file order,
// e.g. BILL v3:  C1C L1C D1C S1C C1W S1W C2W L2W D2W S2W C2L L2L D2L ...
//          index:  0   1   2   3   4   5   6   7   8   9  10  11  12
//
// Each satellite's data line holds its values in that SAME order. So to
// fill an ObsRecord we only need four numbers: the list position of the
// L1 code, L2 code, L1 phase and L2 phase we want to use. For BILL v3:
//   P1 = C1W -> 4,   P2 = C2W -> 6,   L1 = L1C -> 1,   L2 = L2W -> 7
// The body readers (parts 7 and 8) then read value #4, #6, #1 and #7 for
// every satellite at every epoch -- the search happens once, not per line.
//
// Indices stores those four positions; -1 means "not in this file".
// ============================================================
struct Indices {
    int p1 = -1;    // L1 code  -> ObsRecord::pseudorange
    int p2 = -1;    // L2 code  -> ObsRecord::pseudorange_c2
    int l1 = -1;    // L1 phase -> ObsRecord::phase_l1
    int l2 = -1;    // L2 phase -> ObsRecord::phase_l2
};

// ------------------------------------------------------------
// pickObs: find the best available code for ONE field.
//
//   types       the header list (file order), e.g. {"C1C","L1C",...,"C1W",...}
//   candidates  acceptable codes, BEST FIRST,  e.g. {"C1W","C1P","C1Y","C1C"}
//   chosen      OUTPUT: the code that won, or "-" if none was found
//   returns     the winner's position in types, or -1 if none was found
//
// The OUTER loop runs over the candidates (priority order) and the INNER
// loop searches the header. That order is essential. 
// Walk-through for P1 on BILL v3: candidates {"C1W","C1P","C1Y","C1C"}:
//   1st candidate "C1W": search the list ... found at index 4 -> return 4
//   (C1C, which is at index 0, is never reached because C1W was found first)
// If the loops were swapped, the first code IN THE FILE would win, and
// C1C (index 0) would beat C1W (index 4) just because it is listed first.
//
// chosen is a reference (string&), not a pointer: the caller ALWAYS wants
// to know which code won, so the output is mandatory -- unlike the optional
// RinexObsHeader* hdr = nullptr in parseRinexAllEpochs.
// ------------------------------------------------------------
static int pickObs(const vector<string>& types, const vector<string>& candidates,
    string& chosen) {
    for (const auto& c : candidates)                          // best candidate first
        for (size_t i = 0; i < types.size(); i++)             // search the header list for it
            if (types[i] == c) {                              // found:
                chosen = c;                                   //   report which code is used
                return static_cast<int>(i);                   //   and its position (size_t -> int)
            }
    chosen = "-";                                             // no candidate exists in this file
    return -1;
}

// ------------------------------------------------------------
// resolveIndices: choose the code for each of the four fields, by version.
// Also writes the winners into h.code_p1 ... h.code_l2, so callers can
// print them and check that base and rover use the SAME signals.
//
// How to read a RINEX 3 code:  <type><band><attribute>
//   type       C = code (pseudorange), L = carrier phase, D = Doppler, S = signal strength
//   band       1 = L1 (1575.42 MHz), 2 = L2 (1227.60 MHz), 5 = L5
//   attribute  tracking mode / signal component:
//                C = C/A,  P = P (no anti-spoofing),  Y = Y,
//                W = semi-codeless P(Y) (what geodetic receivers record under A/S),
//                L / S / X = L2C (L, M or combined)
//   e.g. C1W = L1 P(Y) code,  L2W = L2 P(Y) phase,  L2L = L2C phase
//
// Priority for each field and why:
//   P1 (L1 code):  P(Y) family first -- C1W, C1P, C1Y -- because the broadcast
//                  clocks and TGD refer to P(Y). Last resort: C1C (C/A code).
//                  C1C works, but P1 and C1 differ by a satellite-specific bias
//                  (P1-C1 DCB, typically up to ~2 ns = tens of cm), which does NOT cancel in double
//                  differences if base and rover pick different codes.
//   P2 (L2 code):  P(Y) family only -- C2W, C2P, C2Y -- the code that matches
//                  the L2 P(Y) phase.
//   L1 (L1 phase): L1C first (what Week 15 used; virtually every receiver
//                  records it), then L1W / L1P.
//   L2 (L2 phase): P(Y) family only -- L2W, L2P, L2Y.
//
// Why L2C codes (C2L/C2S/C2X, L2L/L2S/L2X) are NOT fallbacks:
//   L2C is a different signal with its own tracking loop, so its phase has a
//   different integer ambiguity and different biases than L2 P(Y). Verified on
//   BILL, G06, first epoch: L2L - L2W = 8.011 cycles. If one station gave us
//   L2W and the other L2L, the L2 double differences would no longer be integers.
//   Better to report "no L2" than to silently mix signals.
//
// RINEX 2 caveat:
//   2-character codes carry no tracking information. "L2" is whatever the
//   converter put there -- in EarthScope's bill0070.24o it is the L2C phase
//   (equal to the v3 L2L value), not L2 P(Y). The parser cannot detect this;
//   only a comparison like the one above (or the receiver documentation) can.
//
// Expected results:
//   BILL v3       P1=C1W[4]  P2=C2W[6]  L1=L1C[1]  L2=L2W[7]
//   bill0070.24o  P1=P1[4]   P2=P2[3]   L1=L1[0]   L2=L2[1]
// ------------------------------------------------------------
static Indices resolveIndices(RinexObsHeader& h) {
    Indices ix;                                               // all four start at -1
    const auto& t = h.gps_obs_types;                          // short name for the header list
    if (h.major == 2) {
        ix.p1 = pickObs(t, { "P1", "C1" }, h.code_p1);          // P(Y) code, else C/A code
        ix.p2 = pickObs(t, { "P2" }, h.code_p2);          // P(Y) code only
        ix.l1 = pickObs(t, { "L1" }, h.code_l1);          // the only L1 phase in v2
        ix.l2 = pickObs(t, { "L2" }, h.code_l2);          // the only L2 phase in v2 (tracking unknown!)
    }
    else {                                                  // RINEX 3 (and 4) obs
        ix.p1 = pickObs(t, { "C1W", "C1P", "C1Y", "C1C" }, h.code_p1);   // P(Y) family, else C/A
        ix.p2 = pickObs(t, { "C2W", "C2P", "C2Y" }, h.code_p2);   // P(Y) family only
        ix.l1 = pickObs(t, { "L1C", "L1W", "L1P" }, h.code_l1);   // C/A phase first
        ix.l2 = pickObs(t, { "L2W", "L2P", "L2Y" }, h.code_l2);   // P(Y) family only
    }
    // {"P1", "C1"}: the braces build a temporary vector<string> right in the call.
    return ix;
}

// ============================================================
// Reading values: the "slot"
// ============================================================
//
// In BOTH RINEX 2 and RINEX 3, every observation value is written in the
// same 16-character slot (Fortran format F14.3, I1, I1):
//
//   chars 0-13   the value, 14 wide, 3 decimals  (pseudorange in m, phase in cycles)
//   char  14     LLI  = loss-of-lock indicator   (0-7, or blank)
//   char  15     SSI  = signal-strength indicator (1-9, or blank; not used here)
//
//   Example slot (L1 phase, LLI = 0, SSI = 8):
//     "  130321269.801" + "0" + "8"
//      |<--- value -->|  LLI  SSI
//
// What differs between versions is only HOW the slots are put on the page,
// and the body readers (parts 7 and 8) undo that before calling readSlot:
//
//   RINEX 3: one line per satellite, "G06" + all slots in a row
//            -> data = the line without its first 3 characters
//   RINEX 2: the slots wrap after 5 per line (80 columns), so one satellite
//            with 20 types takes 4 lines
//            -> data = those 4 lines glued together, each padded to 80 chars
//
// Either way, slot i then starts at character 16*i of 'data', so one
// function reads both versions:
//
//   data:  [ slot 0 ][ slot 1 ][ slot 2 ][ slot 3 ] ...
//   index:  0         16        32        48
//
// Missing values: RINEX writes a missing observation as blanks (both
// versions) or as 0.0 (RINEX 2). Blanks -> toDouble returns its fallback 0.0,
// so "value == 0.0" catches both cases. A real pseudorange or phase is never
// exactly 0.000.
//
// LLI bits (we keep the whole digit; callers test the bit they need):
//   bit 0 (value 1): lost lock -> a cycle slip is possible   <- used as  lli & 1
//   bit 1 (value 2): half-cycle ambiguity possible
//   bit 2 (value 4): RINEX 2: observed under anti-spoofing   (NOT a slip; very
//                    common on L2 in RINEX 2 files -- bill0070.24o says so in
//                    its header: "BIT 2 OF LLI FLAGS DATA COLLECTED UNDER A/S")
// ============================================================
struct Slot {
    double value = 0.0;   // the observation (m or cycles)
    int    lli = 0;     // LLI digit, 0 if blank
    bool   ok = false; // true only if a real (non-zero) value was found
};

// Read slot number idx from one satellite's joined data string.
static Slot readSlot(const string& data, int idx) {
    Slot s;                                                   // starts as "missing"
    if (idx < 0) return s;                                    // this observable is not in the file
    const size_t pos = 16 * static_cast<size_t>(idx);         // first character of slot idx
    const double v = toDouble(field(data, pos, 14), 0.0);     // the 14-char value; blank -> 0.0
    if (v == 0.0) return s;                                   // blank or 0.0 = missing observation
    s.value = v;
    s.ok = true;
    const char c = field(data, pos + 14, 1)[0];               // the LLI character (blank if line was cut)
    s.lli = (c >= '0' && c <= '9') ? c - '0' : 0;             // '4' -> 4 ; blank -> 0
    return s;                                                 // ('4' - '0' = 4: digit chars are consecutive)
}

// ============================================================
// Build one ObsRecord from one satellite's data string.
//
//   prn   satellite number (6 for G06)
//   data  that satellite's joined slots (see above)
//   ix    which slot holds P1, P2, L1, L2 (from resolveIndices)
//   r     OUTPUT: the filled record
//   returns false if the satellite has no L1 code -> caller drops it
//
// Rule kept from Week 15: a satellite is only kept if it has an L1 code.
// Everything else is optional and flagged with has_c2 / has_phase_l1 /
// has_phase_l2, so downstream code (slip detection, DD) can check what exists.
// ============================================================
static bool buildRecord(int prn, const string& data, const Indices& ix, ObsRecord& r) {
    const Slot p1 = readSlot(data, ix.p1);                    // L1 code
    if (!p1.ok) return false;                                 // no L1 code -> skip this satellite
    r = ObsRecord{};                                          // reset every field to its default
    r.prn = prn;
    r.pseudorange = p1.value;

    const Slot p2 = readSlot(data, ix.p2);                    // L2 code (optional)
    if (p2.ok) { r.pseudorange_c2 = p2.value; r.has_c2 = true; }

    const Slot l1 = readSlot(data, ix.l1);                    // L1 phase (optional) + its LLI
    if (l1.ok) { r.phase_l1 = l1.value; r.has_phase_l1 = true; r.lli_l1 = l1.lli; }

    const Slot l2 = readSlot(data, ix.l2);                    // L2 phase (optional) + its LLI
    if (l2.ok) { r.phase_l2 = l2.value; r.has_phase_l2 = true; r.lli_l2 = l2.lli; }

    return true;
}

// ============================================================
// RINEX 3 body
// ============================================================
//
// After END OF HEADER, a RINEX 3 file is a sequence of epoch blocks:
//
//   > 2024 01 07 00 00  0.0000000  0 NN        <- epoch line: time, flag, NN satellites
//   G06  115534040.685 ... (all slots) ...     <- NN satellite lines, one per satellite,
//   E26  ...                                       in any system order
//   R21  ...
//   > 2024 01 07 00 00 15.0000000  0 NN        <- next epoch
//   ...
//
// Epoch line (format A1,1X,I4.4,4(1X,I2.2),F11.7,2X,I1,I3,6X,F15.12):
//   index  0       '>'          epoch marker (only epoch lines start with it)
//          2-5     2024         year (4 digits)
//          7-8     01           month
//          10-11   07           day
//          13-14   00           hour
//          16-17   00           minute
//          18-28   0.0000000    seconds, 11 wide, 7 decimals
//          31      0            epoch flag
//          32-34   NN           number of satellites (or of special records, see below)
//          41-55   (optional)   receiver clock offset -- not used here
//
// Epoch flag (column 32 in 1-based counting):
//   0  OK
//   1  power failure since the previous epoch -- the data are still valid, keep them
//   2-5  EVENTS, not data: 2 = antenna starts moving, 3 = new site occupation,
//        4 = header information follows, 5 = external event.
//        Here NN = number of SPECIAL RECORD LINES that follow -> skip exactly NN lines
//   6  cycle-slip records follow: NN lines that LOOK like satellite lines but
//      are receiver slip reports, not observations -> read them and discard
//
// Satellite line:
//   index  0       G            system letter (G, R, E, C, J, S, I)
//          1-2     06           PRN
//          3 ...                the slots (16 chars each), in header-list order
//   -> line.substr(3) is exactly the 'data' string readSlot expects.
//
// Staying in sync: we read EXACTLY NN lines after every epoch line, including
// the non-GPS ones we throw away. If a line ever surprises us, the outer loop
// skips forward until the next line that starts with '>' (self-repairing).
//
//   file      the open stream, positioned just after END OF HEADER
//   ix        slot positions of P1, P2, L1, L2 (from resolveIndices)
//   epochs    OUTPUT: one ObsEpoch is appended per data epoch
//   n_events  OUTPUT: count of skipped event / cycle-slip blocks (for the summary)
// ============================================================
static void readBodyV3(istream& file, const Indices& ix,
    vector<ObsEpoch>& epochs, int& n_events) {
    string line;
    while (getLine(file, line)) {
        if (line.empty() || line[0] != '>') continue;        // not an epoch line: skip until one appears

        const int flag = toInt(field(line, 31, 1));           // epoch flag (0-6)
        const int n = toInt(field(line, 32, 3));           // satellites, or special-record lines

        if (flag >= 2 && flag <= 5) {                         // event block, no observations:
            for (int k = 0; k < n && getLine(file, line); k++) {}   // skip its n lines
            n_events++;
            continue;                                         // go look for the next '>'
        }

        // ---- Epoch time ----
        ObsEpoch ep;
        ep.year = toInt(field(line, 2, 4));
        ep.month = toInt(field(line, 7, 2));
        ep.day = toInt(field(line, 10, 2));
        ep.hour = toInt(field(line, 13, 2));
        ep.minute = toInt(field(line, 16, 2));
        ep.second = toDouble(field(line, 18, 11));
        toGpsTime(ep);                                        // fills gps_week and t_gps (part 2)

        // ---- The n satellite lines of this epoch ----
        for (int k = 0; k < n; k++) {
            if (!getLine(file, line)) break;                  // file ended early (truncated file)
            if (flag == 6) continue;                          // cycle-slip records: consume, don't use
            if (line.size() < 3 || line[0] != 'G') continue;  // not GPS: line consumed, nothing kept
            const int prn = toInt(field(line, 1, 2));         // "G06" -> 6
            ObsRecord r;
            if (buildRecord(prn, line.substr(3), ix, r))      // slots start after "Gnn"
                ep.gps_records.push_back(r);                  // keep only satellites with an L1 code
        }

        if (flag == 6) { n_events++; continue; }              // slip-record block: no epoch to store
        epochs.push_back(ep);                                 // flags 0 and 1: a real data epoch
    }
}

// ============================================================
// RINEX 2 body
// ============================================================
//
// RINEX 2 differs from RINEX 3 in three ways that matter here:
//   1. No '>' marker. An epoch line is recognized only by POSITION: it is the
//      line right after the previous epoch's data (or after END OF HEADER).
//   2. The satellite list sits ON the epoch line (12 per line, more on
//      continuation lines), and the data lines carry NO satellite ID.
//      Data block k belongs to satellite k of the list -- order is everything.
//   3. Each satellite's slots wrap after 5 per line (5 x 16 = 80 columns),
//      so one satellite takes ceil(n_types / 5) lines.
//
// Example, bill0070.24o (20 types -> 4 lines per satellite, 29 satellites):
//
//    24  1  7  0  0  0.0000000  0 29G06G17G24E26G14G11R21G12G30G13E27E19   <- epoch line, sats 1-12
//                                  E21R02S31R09G15E13R22E0...              <- continuation, sats 13-24
//                                  ...                                     <- continuation, sats 25-29
//    115534040.685   90026595.611   21985377.340 ...                       <- G06, data line 1 of 4
//   ...                                                                    <- G06, data lines 2-4
//   ...                                                                    <- G17, 4 lines, and so on
//   -> one epoch = 3 list lines + 29 x 4 = 116 data lines
//
// Epoch line (format 1X,I2.2,4(1X,I2),F11.7,2X,I1,I3,12(A1,I2)):
//   index  1-2     24           year, 2 digits: 80-99 -> 19xx, 00-79 -> 20xx
//          4-5      1           month
//          7-8      7           day
//          10-11    0           hour
//          13-14    0           minute
//          15-25   0.0000000    seconds, 11 wide, 7 decimals
//          28      0            epoch flag (same meanings as RINEX 3, see part 7)
//          29-31   29           number of satellites (or of special records)
//          32-67   G06G17...    up to 12 satellite IDs, 3 chars each: 32 + 3*j
//          68-79   (optional)   receiver clock offset -- not used here
// Continuation line for satellites 13-24, 25-36, ...:
//          0-31    blanks       (format 32X,12(A1,I2))
//          32-67   E21R02...    the next 12 IDs, same columns as on the epoch line
//
// Satellite ID: system letter + PRN. The letter may be blank in a GPS-only
// file (blank = GPS), and the PRN may be space-padded ("G 6" instead of "G06").
//
// Events: flags 2-5 -> n = number of special-record lines to skip (their time
// fields may even be blank). Flag 6 -> a satellite list + data block in the
// normal layout that holds slip reports, not observations -> read, discard.
//
//   file      the open stream, positioned just after END OF HEADER
//   ix        slot positions of P1, P2, L1, L2 (from resolveIndices)
//   n_types   number of observation types in the header (20 for BILL)
//   epochs    OUTPUT: one ObsEpoch is appended per data epoch
//   n_events  OUTPUT: count of skipped event / cycle-slip blocks
// ============================================================
static void readBodyV2(istream& file, const Indices& ix, int n_types,
    vector<ObsEpoch>& epochs, int& n_events) {
    const int lines_per_sat = (n_types + 4) / 5;             // integer ceiling of n_types/5: 20 -> 4, 7 -> 2
    string line;                                              // the EPOCH line (kept: its time is read later)
    string cur;                                               // scratch: continuation and data lines

    while (getLine(file, line)) {
        if (trim(line).empty()) continue;                     // tolerate blank lines between epochs

        const int flag = toInt(field(line, 28, 1));           // epoch flag (0-6)
        const int n = toInt(field(line, 29, 3));           // satellites, or special-record lines

        if (flag >= 2 && flag <= 5) {                         // event block, no observations:
            for (int k = 0; k < n && getLine(file, cur); k++) {}    // skip its n lines
            n_events++;
            continue;
        }

        // ---- 1. Satellite list: 12 IDs per line, continuation lines for more ----
        vector<string> sats;
        cur = line;                                           // IDs 1-12 are on the epoch line itself
        for (int k = 0; k < n; k++) {
            if (k > 0 && k % 12 == 0 && !getLine(file, cur)) return;   // 13th, 25th, ... ID: next line
            sats.push_back(field(cur, 32 + 3 * (k % 12), 3));          // ID j of its line at 32 + 3j
        }

        // ---- 2. Epoch time (from 'line', which still holds the epoch line) ----
        ObsEpoch ep;
        const int yy = toInt(field(line, 1, 2));
        ep.year = (yy < 80) ? 2000 + yy : 1900 + yy;        // RINEX 2 two-digit-year rule
        ep.month = toInt(field(line, 4, 2));
        ep.day = toInt(field(line, 7, 2));
        ep.hour = toInt(field(line, 10, 2));
        ep.minute = toInt(field(line, 13, 2));
        ep.second = toDouble(field(line, 15, 11));
        toGpsTime(ep);                                        // fills gps_week and t_gps (part 2)

        // ---- 3. Data blocks, in the SAME order as the satellite list ----
        for (const string& s : sats) {
            string data;                                      // this satellite's slots, joined
            for (int L = 0; L < lines_per_sat; L++) {         // read ALL its lines, even for R/E/S,
                if (!getLine(file, cur)) break;               //   so the next block starts in the right place
                data += field(cur, 0, 80);                    // pad each line back to 80 (trailing blanks may be cut)
            }                                                 //   -> slot i is at 16*i, as readSlot expects
            if (flag == 6) continue;                          // slip-report block: consume, don't use
            if (s[0] != 'G' && s[0] != ' ') continue;         // keep GPS only (blank letter = GPS)
            const int prn = toInt(s.substr(1, 2));            // "G06" or "G 6" -> 6
            ObsRecord r;
            if (buildRecord(prn, data, ix, r))
                ep.gps_records.push_back(r);                  // keep only satellites with an L1 code
        }

        if (flag == 6) { n_events++; continue; }              // slip-record block: no epoch to store
        epochs.push_back(ep);                                 // flags 0 and 1: a real data epoch
    }
}

// ============================================================
// Public entry point
// ============================================================
//
// The whole parser in one picture:
//
//   parseRinexAllEpochs(filename, &h)
//     |
//     |-- detectRinexVersion()   line 1 only: exists? compressed? CRINEX? RINEX?
//     |                          (prints the reason and stops if not usable)
//     |-- readHeader()           version, marker, XYZ, interval, GPS type list
//     |-- resolveIndices()       which slot holds P1, P2, L1, L2 (once per file)
//     |-- print a summary        version, reader used, chosen codes, warnings
//     |-- readBodyV2() or readBodyV3()   every epoch -> vector<ObsEpoch>
//     |-- copy the header to the caller (only if they passed &h)
//     '-- return the epochs
//
// Note: the parameter is called hdr_out here and hdr in rinex_obs.h. That is
// allowed -- only the TYPES must match between declaration and definition.
// The default "= nullptr" is written in the header only (part 3.1).
//
// Returning a vector by value is cheap in modern C++: the compiler moves it
// to the caller instead of copying 5760 epochs.
// ============================================================
vector<ObsEpoch> parseRinexAllEpochs(const string& filename, RinexObsHeader* hdr_out) {
    vector<ObsEpoch> epochs;                                  // stays empty if anything fails

    if (detectRinexVersion(filename) == 0.0) return epochs;   // prints why (missing, compressed, CRINEX, ...)

    ifstream file(filename);                                  // open again for the real read
    if (!file.is_open()) {
        cerr << "Cannot open: " << filename << endl;
        return epochs;
    }

    // ---- Header ----
    RinexObsHeader h;
    if (!readHeader(file, h)) return epochs;                  // readHeader printed the reason
    if (h.gps_obs_types.empty()) {                            // e.g. a Galileo-only file
        cerr << "WARNING: no GPS observation types in header" << endl;
        return epochs;
    }

    // ---- Which slots to read ----
    const Indices ix = resolveIndices(h);                     // also fills h.code_p1 ... h.code_l2

    // ---- Summary ----
    // Format the version in a local string stream, so cout's own settings
    // (fixed, precision) are NOT changed for the caller's later output.
    ostringstream ver;
    ver << fixed << setprecision(2) << h.version;             // 3.04, not 3.0399999
    cout << "RINEX " << ver.str()
        << " -> v" << (h.major == 2 ? 2 : 3) << " reader"    // 4.xx obs also uses the v3 reader
        << "  | marker: " << (h.marker_name.empty() ? "?" : h.marker_name) << endl;
    cout << "Resolved indices: P1=" << h.code_p1 << "[" << ix.p1 << "]"
        << "  P2=" << h.code_p2 << "[" << ix.p2 << "]"
        << "  L1=" << h.code_l1 << "[" << ix.l1 << "]"
        << "  L2=" << h.code_l2 << "[" << ix.l2 << "]" << endl;
    if (ix.p1 < 0) cerr << "WARNING: no L1 code observable -> no records will be kept" << endl;
    if (ix.p2 < 0) cerr << "WARNING: no L2 P(Y) code -> dual-frequency unavailable" << endl;
    if (ix.l2 < 0) cerr << "WARNING: no L2 P(Y) phase -> LG/MW unavailable" << endl;

    // ---- Body: the only place the two versions take different paths ----
    int n_events = 0;
    if (h.major == 2) readBodyV2(file, ix, (int)h.gps_obs_types.size(), epochs, n_events);
    else              readBodyV3(file, ix, epochs, n_events);

    if (n_events > 0)
        cout << "Skipped " << n_events << " event-flag records (flags 2-6)" << endl;

    if (hdr_out) *hdr_out = h;                                // optional output (see part 3.1)
    return epochs;
}