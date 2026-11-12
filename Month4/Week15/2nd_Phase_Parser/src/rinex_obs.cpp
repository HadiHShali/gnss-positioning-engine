// ============================================================================
//  rinex_obs.cpp
// ----------------------------------------------------------------------------
//  The bodies of the functions declared in rinex_obs.h: read a RINEX 3
//  observation file and fill ObsEpoch / ObsRecord (GPS satellites only).
//
//  Order of this file: small tools first, then the functions that use them.
//
//      Section 1   includes
//      Section 2   time helpers       calendar date -> GPS week, sec of week
//      Section 3   readObsField       read one 16-character observation slot
//      Section 4   findObsIndex       which position holds "L1C"?
//      Section 5   readEpochs         the ONE reader behind both public
//                                     functions (Parts A to D)
//      Section 6   parseRinexAllEpochs and parseRinexObsOneEpoch
//                                     the two public functions, short
// ============================================================================


// ============================================================================
//  SECTION 1 -- INCLUDES
// ----------------------------------------------------------------------------
//  "rinex_obs.h" : comes first. The compiler then checks that the header
//                  works on its own, and that the bodies below match its
//                  declarations.
//  <cctype>      : std::isdigit, to test the LLI character.
//  <fstream>     : std::ifstream, to open and read the file.
//  <iostream>    : std::cout and std::cerr, for messages.
//  <string> and <vector> arrive through rinex_obs.h.
//  using namespace std : acceptable in a .cpp, because it affects this file
//                  only. NEVER write it in a header: it would leak into
//                  every file that includes that header.
// ============================================================================
#include "rinex_obs.h"
#include <cctype>
#include <fstream>
#include <iostream>

using namespace std;


// ============================================================================
//  SECTION 2 -- TIME HELPERS : CALENDAR DATE -> GPS WEEK + SECONDS OF WEEK
// ----------------------------------------------------------------------------
//  Why functions : pure calculations. A date goes in, a number comes out.
//  Why "static"  : makes a function private to this .cpp. No other file can
//                  call it or clash with its name. A helper that is not in
//                  the header should be static.
//  Time scale    : the epoch times in a GPS RINEX file are already in GPS
//                  time, so no leap seconds are applied here.
//  GPS week      : week 0 began on Sunday 1980-01-06. Seconds of week run
//                  from 0 (Sunday 00:00:00) up to 604800.
// ============================================================================

// GPS week number = whole weeks elapsed since 1980-01-06.
static int computeGpsWeek(int year, int month, int day) {
    // Julian Day Number of the date (standard integer formula).
    int a   = (14 - month) / 12;            // 1 for January/February, else 0
    int y   = year + 4800 - a;              // shifted year
    int m   = month + 12 * a - 3;           // March = 0 ... February = 11
    int jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400
            - 32045;

    const int JDN_GPS_START = 2444245;      // Julian Day Number of 1980-01-06
    return (jdn - JDN_GPS_START) / 7;       // integer division = whole weeks
}

// Seconds since the start of the GPS week (Sunday 00:00:00).
static double computeGpsSecOfWeek(int year, int month, int day,
                                  int hour, int min, double sec) {
    // Zeller's congruence gives the day of the week.
    // January and February count as months 13 and 14 of the previous year.
    int y = year, m = month;
    if (m < 3) { m += 12; y -= 1; }
    int K = y % 100;                        // year within the century
    int J = y / 100;                        // century

    // h: 0 = Saturday, 1 = Sunday, 2 = Monday, ..., 6 = Friday
    int h = (day + (13 * (m + 1)) / 5 + K + K / 4 + J / 4 + 5 * J) % 7;

    // GPS day of week: 0 = Sunday, 1 = Monday, ..., 6 = Saturday
    int gps_dow = (h + 6) % 7;

    return gps_dow * 86400.0 + hour * 3600.0 + min * 60.0 + sec;
}


// ============================================================================
//  SECTION 3 -- readObsField : READ ONE 16-CHARACTER OBSERVATION SLOT
// ----------------------------------------------------------------------------
//  Why it exists : the same steps (cut 14 characters, check for blank,
//                  convert, read the LLI) were written three times, for C2W,
//                  L1C and L2W. Repeated code is a function waiting to be
//                  written: a fix made in one copy is easily forgotten in
//                  the others.
//  Why a struct  : one slot gives three values that travel together
//                  (value, LLI, present), so the function returns them as
//                  one ObsField.
//  namespace { } : everything inside is private to this .cpp. It is the
//                  same idea as "static", and it also works for a struct.
//  Slot layout   : value = columns col .. col+13   (F14.3)
//                  LLI   = column  col+14
//                  SSI   = column  col+15          (not used here)
//  "present" stays false when:
//                  - the observable is not in this file      (col < 0)
//                  - the line ends before the slot  (trailing blanks are
//                    often cut off)
//                  - the field is blank
//                  - the value is 0.0   (RINEX: missing = blank OR 0.0)
// ============================================================================
namespace {

struct ObsField {
    double value   = 0.0;                   // the number, in the file's units
    int    lli     = 0;                     // loss-of-lock digit, 0 if blank
    bool   present = false;                 // true only for a real observation
};

ObsField readObsField(const string& line, int col) {
    ObsField f;                                         // starts "not present"

    if (col < 0) return f;                              // not in this file
    if ((int)line.size() < col + 14) return f;          // line too short

    try {
        f.value = stod(line.substr(col, 14));           // text -> number
    } catch (...) {
        return f;                                       // blank or damaged
    }
    if (f.value == 0.0) return f;                       // 0.0 means missing

    f.present = true;                                   // a real observation

    int lli_col = col + 14;                             // LLI follows the value
    if ((int)line.size() > lli_col &&
        isdigit(static_cast<unsigned char>(line[lli_col]))) {
        f.lli = line[lli_col] - '0';                    // character -> number
    }
    return f;
}

}  // namespace


// ============================================================================
//  SECTION 4 -- findObsIndex : WHICH POSITION HOLDS A GIVEN OBSERVABLE CODE?
// ----------------------------------------------------------------------------
//  This is the body of the function declared in rinex_obs.h. The signature
//  here must match that declaration exactly.
//  Not static    : it is in the header, so other files may call it.
//  How           : walk through the list and return the position of the
//                  first match, or -1 if the code is not there.
//  size_t vs int : obs_types.size() is a size_t (unsigned), which cannot
//                  hold -1. The function returns int so that it can, and
//                  static_cast<int> makes the conversion explicit.
// ============================================================================
int findObsIndex(const vector<string>& obs_types, const string& code) {
    for (size_t i = 0; i < obs_types.size(); i++) {
        if (obs_types[i] == code) return static_cast<int>(i);   // found
    }
    return -1;                                                  // not found
}


// ============================================================================
//  SECTION 5 -- readEpochs : THE ONE READER BEHIND BOTH PUBLIC FUNCTIONS
// ----------------------------------------------------------------------------
//  Why one shared function : the two public functions used to be written
//      separately, and they drifted apart. One read C1C at a fixed column,
//      the other read C1W; one used seconds of week, the other seconds of
//      day. With a single body they cannot disagree.
//  Why "static" : it is a private helper, not declared in the header.
//  Inputs :
//      filename    the RINEX 3 observation file
//      min_t_gps   skip epochs earlier than this (seconds of GPS week);
//                  0 means "start from the first epoch"
//      max_epochs  stop after storing this many epochs; 0 means "no limit"
//  Output : the epochs that were read; empty if the file could not be used.
//  The function does four jobs in order. Each has its own box below:
//      Part A   open the file, read the GPS observable list from the header
//      Part B   turn the wanted codes into column numbers
//      Part C   read one epoch line (flag, count, time)
//      Part D   read that epoch's satellite lines, then store the epoch
// ============================================================================
static vector<ObsEpoch> readEpochs(const string& filename,
                                   double min_t_gps, int max_epochs) {
    vector<ObsEpoch> epochs;                // the result; empty on failure

    // ------------------------------------------------------------------------
    //  PART A -- OPEN THE FILE AND READ THE GPS OBSERVABLE LIST
    // ------------------------------------------------------------------------
    //  The header record "SYS / # / OBS TYPES" lists, for each constellation,
    //  the observables in the order they appear in every satellite line:
    //      column 0        system letter ('G' = GPS)
    //      columns 3-5     how many observables
    //      from column 7   the 3-character codes, one every 4 columns,
    //                      at most 13 per line; the rest continue on the
    //                      next line, again starting at column 7.
    //  Version check : this reader understands RINEX 3 (and 4, same layout).
    //                  A RINEX 2 file has a different layout, so it is
    //                  refused with a message instead of giving zero epochs.
    //  try/catch     : a damaged header stops this file with a message; it
    //                  does not end the program, so a multi-day run can go
    //                  on to the next file.
    // ------------------------------------------------------------------------
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Cannot open OBS file: " << filename << endl;
        return epochs;
    }

    string line;
    vector<string> gps_obs_types;           // for example {"C1C","L1C",...}

    try {
        while (getline(file, line)) {
            if (line.find("END OF HEADER") != string::npos) break;  // done

            // ---- first header line: the format version, columns 0-8 ----
            if (line.find("RINEX VERSION / TYPE") != string::npos) {
                double version = stod(line.substr(0, 9));   // e.g. 3.04
                if (version < 3.0) {
                    cerr << "RINEX " << version << " is not supported "
                         << "(RINEX 3 or later needed): " << filename << endl;
                    return epochs;
                }
            }

            // ---- the GPS observable list ----
            bool is_obs_types =
                line.find("SYS / # / OBS TYPES") != string::npos;
            if (is_obs_types && line[0] == 'G') {
                int n_obs = stoi(line.substr(3, 3));    // how many observables
                gps_obs_types.clear();

                int    col      = 7;                // first code starts here
                int    parsed   = 0;                // codes read so far
                string cur_line = line;             // the line being read

                while (parsed < n_obs) {
                    if (col + 3 > 60) {             // past the 13th code:
                        getline(file, cur_line);    //   take the next line
                        col = 7;                    //   codes start at 7 again
                    }
                    string code = cur_line.substr(col, 3);      // e.g. "L1C"
                    code.erase(code.find_last_not_of(' ') + 1); // trim right
                    gps_obs_types.push_back(code);
                    col += 4;                       // 3 characters + 1 blank
                    parsed++;
                }
            }
        }
    } catch (...) {
        cerr << "Damaged header in: " << filename << endl;
        return epochs;
    }

    // ------------------------------------------------------------------------
    //  PART B -- TURN THE WANTED CODES INTO COLUMN NUMBERS
    // ------------------------------------------------------------------------
    //  findObsIndex gives the slot number of a code (or -1). obsColumn turns
    //  a slot number into a column. A missing code keeps column -1, and
    //  readObsField treats -1 as "not present".
    //  C1 fallback : if the file has no C1W, use C1C instead and say so.
    //  With NO L1 pseudorange at all nothing can be stored, so stop here
    //  with a message instead of returning an empty result silently.
    // ------------------------------------------------------------------------
    string c1_code = "C1W";
    int idx_C1 = findObsIndex(gps_obs_types, c1_code);
    if (idx_C1 < 0) {                               // no C1W: try C1C
        c1_code = "C1C";
        idx_C1  = findObsIndex(gps_obs_types, c1_code);
    }
    int idx_C2 = findObsIndex(gps_obs_types, "C2W");
    int idx_L1 = findObsIndex(gps_obs_types, "L1C");
    int idx_L2 = findObsIndex(gps_obs_types, "L2W");

    if (idx_C1 < 0) {
        cerr << "No C1W or C1C observable in: " << filename << endl;
        return epochs;
    }

    int col_C1 = obsColumn(idx_C1);                       // valid, checked above
    int col_C2 = (idx_C2 >= 0) ? obsColumn(idx_C2) : -1;  // -1 = not in file
    int col_L1 = (idx_L1 >= 0) ? obsColumn(idx_L1) : -1;
    int col_L2 = (idx_L2 >= 0) ? obsColumn(idx_L2) : -1;

    cout << "Resolved slots: " << c1_code << "=" << idx_C1
         << " C2W=" << idx_C2 << " L1C=" << idx_L1 << " L2W=" << idx_L2 << endl;
    if (col_L1 < 0 || col_L2 < 0) {
        cerr << "Warning: L1C or L2W is missing, carrier phase unavailable"
             << endl;
    }

    // ------------------------------------------------------------------------
    //  PART C -- READ ONE EPOCH LINE
    // ------------------------------------------------------------------------
    //  Every epoch starts with a line beginning with '>'. Layout, columns
    //  counted from 0:
    //      2-5     year          7-8    month       10-11  day
    //      13-14   hour          16-17  minute      18-28  seconds (F11.7)
    //      31      epoch flag    32-34  number of lines that follow
    //  Epoch flag : 0 = normal
    //               1 = power failure since the previous epoch (expect a
    //                   slip on every satellite)
    //               2 or more = special event. The lines that follow are
    //                   NOT satellites, so they are skipped.
    //  try/catch  : stoi, stod and substr throw on a damaged or short line.
    //               Catching here skips that epoch instead of ending the
    //               program.
    // ------------------------------------------------------------------------
    while (getline(file, line)) {
        if (line.empty() || line[0] != '>') continue;   // not an epoch line

        int flag = 0, n_sats = 0;
        try {
            flag   = stoi(line.substr(31, 1));          // epoch flag
            n_sats = stoi(line.substr(32, 3));          // lines that follow
        } catch (...) {
            continue;                                   // damaged epoch line
        }

        if (flag > 1) {                                 // special event:
            for (int i = 0; i < n_sats; i++) {
                getline(file, line);                    // skip its lines
            }
            continue;
        }

        ObsEpoch epoch;
        epoch.epoch_flag = flag;
        try {
            epoch.year   = stoi(line.substr(2, 4));
            epoch.month  = stoi(line.substr(7, 2));
            epoch.day    = stoi(line.substr(10, 2));
            epoch.hour   = stoi(line.substr(13, 2));
            epoch.minute = stoi(line.substr(16, 2));
            epoch.second = stod(line.substr(18, 11));
        } catch (...) {
            continue;                                   // damaged date or time
        }

        // The same instant as GPS week + seconds of week, as the header says.
        epoch.gps_week = computeGpsWeek(epoch.year, epoch.month, epoch.day);
        epoch.t_gps    = computeGpsSecOfWeek(epoch.year, epoch.month, epoch.day,
                                             epoch.hour, epoch.minute,
                                             epoch.second);

        // Too early: skip this epoch. Its satellite lines do not start
        // with '>', so the test at the top of the loop passes over them.
        if (epoch.t_gps < min_t_gps) continue;

        // --------------------------------------------------------------------
        //  PART D -- READ THIS EPOCH'S SATELLITE LINES, THEN STORE THE EPOCH
        // --------------------------------------------------------------------
        //  n_sats lines follow, one per satellite of ANY constellation.
        //  Keep the lines that start with 'G' (GPS) and skip the rest.
        //  Each wanted observable is read by readObsField (Section 3), so
        //  blank fields, 0.0 values and short lines are handled in one place.
        //  A satellite is stored only if it has an L1 pseudorange. C2W and
        //  the two phases are optional and carry their own has_ flags.
        // --------------------------------------------------------------------
        for (int i = 0; i < n_sats; i++) {
            if (!getline(file, line)) break;            // file ended early
            if (line.size() < 3 || line[0] != 'G') continue;   // not GPS

            ObsRecord rec;
            try {
                rec.prn = stoi(line.substr(1, 2));      // "G06" -> 6
            } catch (...) {
                continue;                               // damaged satellite ID
            }

            ObsField c1 = readObsField(line, col_C1);   // L1 pseudorange, m
            if (!c1.present) continue;                  // required: skip sat
            rec.pseudorange = c1.value;

            ObsField c2 = readObsField(line, col_C2);   // L2 pseudorange, m
            rec.pseudorange_c2 = c2.value;
            rec.has_c2         = c2.present;

            ObsField l1 = readObsField(line, col_L1);   // L1 phase, cycles
            rec.phase_l1     = l1.value;
            rec.has_phase_l1 = l1.present;
            rec.lli_l1       = l1.lli;

            ObsField l2 = readObsField(line, col_L2);   // L2 phase, cycles
            rec.phase_l2     = l2.value;
            rec.has_phase_l2 = l2.present;
            rec.lli_l2       = l2.lli;

            epoch.gps_records.push_back(rec);           // keep this satellite
        }

        if (!epoch.gps_records.empty()) {
            epochs.push_back(epoch);                    // keep non-empty epochs
        }

        if (max_epochs > 0 && (int)epochs.size() >= max_epochs) {
            break;                                      // enough: stop reading
        }
    }

    return epochs;
}


// ============================================================================
//  SECTION 6 -- THE TWO PUBLIC FUNCTIONS
// ----------------------------------------------------------------------------
//  These are the bodies of the two functions declared in rinex_obs.h.
//  Each one only chooses the settings and lets readEpochs (Section 5) do
//  the work. Short public functions over one private worker is a pattern
//  worth reusing: the reading rules exist in exactly one place.
//
//  parseRinexAllEpochs   : no minimum time, no limit -> every epoch.
//  parseRinexObsOneEpoch : first epoch at or after target_t_gps, then stop.
//                          A target of 0 gives the first epoch in the file.
//                          Returns an empty ObsEpoch (no satellites) if no
//                          such epoch exists, so check gps_records.empty().
//                          target_t_gps is seconds of the GPS week and has
//                          no week number, so it is only unambiguous inside
//                          one daily file. For several days, read each file
//                          with parseRinexAllEpochs and pick epochs by
//                          gps_week + t_gps.
// ============================================================================
vector<ObsEpoch> parseRinexAllEpochs(const string& filename) {
    return readEpochs(filename, 0.0, 0);            // from the start, no limit
}

ObsEpoch parseRinexObsOneEpoch(const string& filename, double target_t_gps) {
    vector<ObsEpoch> found = readEpochs(filename, target_t_gps, 1);  // just one

    if (found.empty()) {
        cerr << "No epoch found" << endl;
        return ObsEpoch{};                          // empty epoch = failure
    }
    return found[0];                                // the single epoch read
}
