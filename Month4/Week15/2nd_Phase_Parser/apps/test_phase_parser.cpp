// ============================================================================
//  test_phase_parser.cpp
// ----------------------------------------------------------------------------
//  Day 2 test program: follow ONE satellite through the data, look for cycle
//  slips with the geometry-free (LG) combination, and print a report that
//  can be read without knowing the code.
//
//  How to run it:
//      test_phase_parser                       the default file (Section 2)
//      test_phase_parser day1.rnx              one file
//      test_phase_parser day1.rnx day2.rnx     several days, in time order
//
//      Section 1   includes
//      Section 2   settings      the numbers you may want to change
//      Section 3   ArcInfo, Event   what the program remembers for the report
//      Section 4   text helpers  times and durations as readable text
//      Section 5   main          Part A  read the file(s)
//                                Part B  pick the satellite, form LG
//                                Part C  new arc, slip, or normal step?
//                                Part D  print the report
// ============================================================================


// ============================================================================
//  SECTION 1 -- INCLUDES
// ----------------------------------------------------------------------------
//  "rinex_obs.h"  : ObsEpoch, ObsRecord, parseRinexAllEpochs, secondsBetween.
//  "cycle_slip.h" : geometryFree, likelySlip, slipThreshold.
//  <cmath>        : std::abs on a double.
//  <iomanip>      : setprecision, setw, setfill, left, right: tidy columns.
//  <iostream>     : cout and cerr.
//  <sstream>      : ostringstream, to build a text string piece by piece.
//  <string>, <vector> : used directly in this file, so included directly.
//  This file needs BOTH headers. cycle_slip.h does not include rinex_obs.h.
// ============================================================================
#include "rinex_obs.h"
#include "cycle_slip.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace std;


// ============================================================================
//  SECTION 2 -- SETTINGS
// ----------------------------------------------------------------------------
//  Why constants at the top : these are the numbers you will want to change.
//                  Named once here, they are easy to find, and the code
//                  below reads as words instead of unexplained numbers.
//  TARGET_PRN    : the one satellite this test follows.
//  MAX_GAP_S     : a time step longer than this starts a NEW ARC. The
//                  receiver has re-acquired the signal, so the ambiguity is
//                  new and comparing LG across the gap means nothing.
//                  60 s is the value used in Navipedia's example.
//  STRONG_FACTOR : a jump of at least this many times the limit is called
//                  STRONG even without a receiver flag. With 2.0 an L1-only
//                  slip (19.0 cm) and an L2-only slip (24.4 cm) are strong,
//                  and only the smallest slip (5.4 cm) lands in WEAK, which
//                  is also where false alarms are.
//  PREVIEW_ROWS  : how many of the first LG values to print, as a quick
//                  check that the numbers look sensible.
//  DEFAULT_FILE  : used when no file name is given on the command line.
// ============================================================================
const int    TARGET_PRN    = 6;             // follow G06
const double MAX_GAP_S     = 60.0;          // seconds
const double STRONG_FACTOR = 2.0;           // times the limit
const int    PREVIEW_ROWS  = 5;             // rows of LG to show
const string DEFAULT_FILE  =
    "../../data/BILL00USA_R_20240070000_01D_15S_MO.rnx";


// ============================================================================
//  SECTION 3 -- ArcInfo AND Event : WHAT IS REMEMBERED FOR THE REPORT
// ----------------------------------------------------------------------------
//  Why structs : each one is a group of plain values that travel together,
//                the same reason as ObsRecord.
//  Why they exist : the loop only COLLECTS what it finds. The report is
//                printed afterwards, from these two lists. Keeping "find"
//                and "print" apart lets the report be laid out as tables.
//  ArcInfo     : one unbroken stretch of tracking.
//  Event       : something worth reporting INSIDE an arc: LG jumped, or the
//                receiver raised its loss-of-lock flag, or both.
//  "index" members are positions in the list of epochs, so the time can be
//  looked up when printing.
// ============================================================================
struct ArcInfo {
    int    first    = -1;           // index of the arc's first usable epoch
    int    last     = -1;           // index of its last usable epoch
    int    n_epochs = 0;            // usable epochs in the arc
    int    n_slips  = 0;            // LG jumps found inside the arc
    double gap_s    = 0.0;          // time since the previous arc ended
    bool   power_failure = false;   // started by an epoch flag of 1
};

struct Event {
    int    index   = -1;            // which epoch
    int    arc     = 0;             // arc number, counted from 1
    double jump_m  = 0.0;           // LG change from the previous epoch, m
    double dt_s    = 0.0;           // time step of that comparison, s
    double limit_m = 0.0;           // the limit used for that time step, m
    bool   lg_jump = false;         // true if |jump| was over the limit
    bool   flag_l1 = false;         // receiver flag (LLI bit 0) on L1
    bool   flag_l2 = false;         // receiver flag (LLI bit 0) on L2
};


// ============================================================================
//  SECTION 4 -- TEXT HELPERS : TIMES AND DURATIONS AS READABLE TEXT
// ----------------------------------------------------------------------------
//  Why functions : the same formatting is needed in several places. Written
//                  once, it cannot differ between them.
//  Why "static"  : private to this file, like the helpers in rinex_obs.cpp.
//  How           : an ostringstream works like cout, but collects the text
//                  in a string instead of printing it.
//                  setfill('0') + setw(2) prints 5 as "05".
//  dateTime      : "2024-01-07 03:15:30"   full form, for the overview
//  shortTime     : "01-07 03:15:30"        month-day and time, for the
//                  tables, so that every row fits in an 80-column window
//  duration      : "4 h 34 min", "1 min 45 s" or "30 s", whichever fits.
//  flagText      : which frequency carried the receiver's flag.
// ============================================================================
static string dateTime(const ObsEpoch& e) {
    ostringstream s;
    s << setfill('0')
      << setw(4) << e.year   << "-" << setw(2) << e.month  << "-"
      << setw(2) << e.day    << " " << setw(2) << e.hour   << ":"
      << setw(2) << e.minute << ":"
      << setw(2) << static_cast<int>(e.second + 0.5);   // whole seconds
    return s.str();                                     // stream -> string
}

static string shortTime(const ObsEpoch& e) {
    return dateTime(e).substr(5);                       // drop "2024-"
}

static string duration(double seconds) {
    long total = static_cast<long>(seconds + 0.5);      // whole seconds
    long h = total / 3600;                              // full hours
    long m = (total % 3600) / 60;                       // remaining minutes
    long s = total % 60;                                // remaining seconds

    ostringstream out;
    if (h > 0)      out << h << " h " << m << " min";
    else if (m > 0) out << m << " min " << s << " s";
    else            out << s << " s";
    return out.str();
}

static string flagText(bool l1, bool l2) {
    if (l1 && l2) return "L1+L2";
    if (l1)       return "L1";
    if (l2)       return "L2";
    return "none";
}


// ============================================================================
//  SECTION 5 -- main : THE PROGRAM
// ----------------------------------------------------------------------------
//  main only organizes the work. Every calculation is done by a function
//  from one of the two headers:
//      reading        parseRinexAllEpochs     (rinex_obs.h)
//      time step      secondsBetween          (rinex_obs.h)
//      LG             geometryFree            (cycle_slip.h)
//      limit          slipThreshold           (cycle_slip.h)
//      decision       likelySlip              (cycle_slip.h)
//  argc / argv : what was typed on the command line. argc is the number of
//                words, argv[0] is the program name, argv[1] onward are the
//                file names.
//  Return value : 0 = the test ran, 1 = nothing could be read.
// ============================================================================
int main(int argc, char* argv[]) {
    cout << "=== PHASE PARSER + CYCLE SLIP TEST ===" << endl;
    cout << fixed << setprecision(3);       // every number: 3 decimals

    // ------------------------------------------------------------------------
    //  PART A -- READ THE FILE(S)
    // ------------------------------------------------------------------------
    //  Each file is read on its own and its epochs are appended to one long
    //  list. Give the files in time order. Because secondsBetween uses the
    //  GPS week as well, the list can cross midnight, a week or a year.
    //  If nothing at all was read, stop with return code 1. Otherwise the
    //  report would say "0 slips", which looks like a clean result.
    // ------------------------------------------------------------------------
    vector<string> filenames;
    for (int i = 1; i < argc; i++) {
        filenames.push_back(argv[i]);               // from the command line
    }
    if (filenames.empty()) {
        filenames.push_back(DEFAULT_FILE);          // none given: use default
    }

    vector<ObsEpoch> epochs;                        // all days, in order
    for (const string& name : filenames) {
        vector<ObsEpoch> one_file = parseRinexAllEpochs(name);
        cout << "Parsed " << one_file.size() << " epochs from " << name << endl;
        epochs.insert(epochs.end(),                 // append this file
                      one_file.begin(), one_file.end());
    }

    if (epochs.empty()) {
        cerr << "No epochs were read. Nothing to test." << endl;
        return 1;
    }

    // ------------------------------------------------------------------------
    //  PART B -- PICK THE SATELLITE AND FORM LG
    // ------------------------------------------------------------------------
    //  The detector compares each usable epoch with the PREVIOUS usable one,
    //  so two things must be remembered between epochs:
    //      lg_prev     the previous LG value, in meters
    //      prev_index  where that epoch is in the list; -1 = none yet
    //  An epoch is usable when the satellite is present AND both phases were
    //  read. geometryFree must never be called with a missing phase.
    //  The loop uses an index k (not "for each") because positions have to
    //  be remembered.
    //  What is found goes into two lists, arcs and events. Nothing about
    //  them is printed here; Part D prints the report.
    // ------------------------------------------------------------------------
    double lg_prev    = 0.0;                // previous LG, meters
    int    prev_index = -1;                 // previous usable epoch; -1 = none

    vector<ArcInfo> arcs;                   // every arc found
    vector<Event>   events;                 // every jump or flag inside an arc
    int    n_used          = 0;             // usable epochs of this satellite
    int    n_flag_at_start = 0;             // receiver flags at an arc start
    double sampling_s      = 0.0;           // smallest time step seen

    cout << endl << "First LG values (they should change slowly):" << endl;

    for (size_t k = 0; k < epochs.size(); k++) {
        const ObsEpoch& ep = epochs[k];             // the current epoch

        for (const ObsRecord& rec : ep.gps_records) {
            if (rec.prn != TARGET_PRN) continue;                    // other sat
            if (!rec.has_phase_l1 || !rec.has_phase_l2) continue;   // need both

            double lg = geometryFree(rec.phase_l1, rec.phase_l2);   // meters
            n_used++;

            if (n_used <= PREVIEW_ROWS) {           // quick look at the numbers
                cout << "  " << shortTime(ep)
                     << "  L1 " << rec.phase_l1 << " cyc"
                     << "  L2 " << rec.phase_l2 << " cyc"
                     << "  LG " << lg << " m" << endl;
            }

            // ----------------------------------------------------------------
            //  PART C -- NEW ARC, SLIP, OR NORMAL STEP?
            // ----------------------------------------------------------------
            //  Receiver flag : LLI bit 0, tested on EACH frequency.
            //  New arc : there is nothing valid to compare with when
            //              - this is the first usable epoch, or
            //              - the time step is longer than MAX_GAP_S, or
            //              - the epoch flag is 1 (power failure: the receiver
            //                lost every satellite).
            //            A new arc is NOT a slip. It starts a new ArcInfo.
            //  Inside an arc : compare LG with the previous epoch, using the
            //            limit for this time step. A short dropout (30 s,
            //            45 s) is still compared, with its larger limit.
            //            If LG jumped OR the receiver raised a flag, the
            //            epoch is stored as an Event.
            //  Always  : remember this epoch for the next comparison. After a
            //            slip that restarts from the new level.
            // ----------------------------------------------------------------
            bool flag_l1 = (rec.lli_l1 & 1) != 0;   // receiver flag on L1
            bool flag_l2 = (rec.lli_l2 & 1) != 0;   // receiver flag on L2

            bool   new_arc = (prev_index < 0);      // first usable epoch
            double dt      = 0.0;                   // time step, seconds
            if (!new_arc) {
                dt = secondsBetween(epochs[prev_index], ep);
                if (dt > MAX_GAP_S || ep.epoch_flag == 1) new_arc = true;
            }

            if (new_arc) {
                ArcInfo arc;                        // start a new arc
                arc.first         = static_cast<int>(k);
                arc.gap_s         = dt;             // 0 for the very first arc
                arc.power_failure = (ep.epoch_flag == 1);
                arcs.push_back(arc);
                if (flag_l1 || flag_l2) n_flag_at_start++;
            } else {
                double limit   = slipThreshold(dt);             // meters
                bool   lg_jump = likelySlip(lg_prev, lg, limit);

                if (lg_jump || flag_l1 || flag_l2) {        // worth reporting
                    Event ev;
                    ev.index   = static_cast<int>(k);
                    ev.arc     = static_cast<int>(arcs.size());
                    ev.jump_m  = lg - lg_prev;
                    ev.dt_s    = dt;
                    ev.limit_m = limit;
                    ev.lg_jump = lg_jump;
                    ev.flag_l1 = flag_l1;
                    ev.flag_l2 = flag_l2;
                    events.push_back(ev);
                    if (lg_jump) arcs.back().n_slips++;
                }

                if (sampling_s == 0.0 || dt < sampling_s) sampling_s = dt;
            }

            arcs.back().last = static_cast<int>(k);     // the arc grows
            arcs.back().n_epochs++;

            lg_prev    = lg;                        // remember for next time
            prev_index = static_cast<int>(k);
        }
    }

    if (arcs.empty()) {
        cout << "  (none)" << endl << endl
             << "G" << setfill('0') << setw(2) << TARGET_PRN << setfill(' ')
             << " has no epoch with both phases in this data." << endl;
        return 0;
    }

    // ------------------------------------------------------------------------
    //  PART D -- PRINT THE REPORT
    // ------------------------------------------------------------------------
    //  Four blocks, each answering one question:
    //      1. Overview     how long was the satellite tracked?
    //      2. Arcs         when was it tracked, and how did each arc start?
    //      3. Events       what happened inside the arcs, and how sure is it?
    //      4. Bottom line  how many slips, and do LG and the receiver agree?
    //  Jumps and limits are printed in CENTIMETERS, because the slips of
    //  interest are a few centimeters to a few tens of centimeters.
    //  Verdict for an LG jump:
    //      STRONG  the receiver flagged it too, or the jump is at least
    //              STRONG_FACTOR times the limit
    //      WEAK    just over the limit and no receiver flag
    //  FLAG ONLY means the receiver raised a flag but LG did not jump.
    //  setw sets the width of the NEXT item only; left / right say on which
    //  side of that width the text sits.
    // ------------------------------------------------------------------------
    const string line(78, '=');                     // a row of 78 '=' signs

    double tracked_s = 0.0;                         // total time in arcs
    for (const ArcInfo& a : arcs) {
        tracked_s += secondsBetween(epochs[a.first], epochs[a.last]);
    }

    int n_slips = 0, n_strong = 0, n_weak = 0;      // LG jumps
    int n_flag_at_jump = 0, n_flag_only = 0;        // receiver flags
    for (const Event& ev : events) {
        bool flagged = ev.flag_l1 || ev.flag_l2;
        bool big     = abs(ev.jump_m) >= STRONG_FACTOR * ev.limit_m;
        if (ev.lg_jump) {
            n_slips++;
            if (flagged || big) n_strong++; else n_weak++;
            if (flagged) n_flag_at_jump++;
        } else {
            n_flag_only++;
        }
    }

    // ---- block 1: overview ----
    cout << endl << line << endl
         << " CYCLE-SLIP REPORT   satellite G"
         << setfill('0') << setw(2) << TARGET_PRN << setfill(' ') << endl
         << line << endl;
    cout << " Data from : " << dateTime(epochs.front()) << endl;
    cout << " Data to   : " << dateTime(epochs.back()) << endl;
    cout << " Sampling  : " << duration(sampling_s) << endl;
    cout << " Tracked   : " << duration(tracked_s) << " in " << arcs.size()
         << " arcs (" << n_used << " epochs)" << endl;

    // ---- block 2: arcs ----
    cout << endl << " ARCS   (an arc = one unbroken stretch of tracking)"
         << endl;
    cout << "  " << left << setw(3) << "#" << setw(15) << "start"
         << setw(15) << "end" << setw(12) << "length"
         << right << setw(6) << "epochs" << setw(6) << "slips"
         << "  " << left << "started by" << endl;

    for (size_t i = 0; i < arcs.size(); i++) {
        const ArcInfo& a = arcs[i];
        string how = "first appearance";
        if (a.power_failure) how = "power failure flag";
        else if (i > 0)      how = "gap of " + duration(a.gap_s);

        cout << "  " << left << setw(3) << (i + 1)
             << setw(15) << shortTime(epochs[a.first])
             << setw(15) << shortTime(epochs[a.last])
             << setw(12) << duration(secondsBetween(epochs[a.first],
                                                    epochs[a.last]))
             << right << setw(6) << a.n_epochs << setw(6) << a.n_slips
             << "  " << left << how << endl;
    }

    // ---- block 3: events ----
    cout << endl << " EVENTS INSIDE THE ARCS" << endl;
    if (events.empty()) {
        cout << "  none" << endl;
    } else {
        cout << "  " << left << setw(15) << "time" << setw(4) << "arc"
             << right << setw(11) << "LG jump" << setw(9) << "limit"
             << setw(9) << "x limit" << setw(7) << "step"
             << "  " << left << setw(7) << "flag" << "verdict" << endl;

        for (const Event& ev : events) {
            bool   flagged = ev.flag_l1 || ev.flag_l2;
            double ratio   = abs(ev.jump_m) / ev.limit_m;
            string verdict = "FLAG ONLY";
            if (ev.lg_jump) {
                bool big = abs(ev.jump_m) >= STRONG_FACTOR * ev.limit_m;
                verdict  = (flagged || big) ? "STRONG" : "WEAK";
            }

            cout << "  " << left << setw(15) << shortTime(epochs[ev.index])
                 << setw(4) << ev.arc << right << setprecision(1)
                 << showpos << setw(8) << ev.jump_m * 100.0 << noshowpos
                 << " cm" << setw(6) << ev.limit_m * 100.0 << " cm"
                 << setw(9) << ratio
                 << setprecision(0) << setw(5) << ev.dt_s << " s"
                 << "  " << left << setw(7)
                 << flagText(ev.flag_l1, ev.flag_l2) << verdict << endl;
        }
    }

    // ---- block 4: bottom line ----
    cout << endl << " BOTTOM LINE" << endl;
    cout << "  Slips found by the LG test   : " << n_slips << endl;
    cout << "      strong                   : " << n_strong
         << "   (receiver flag, or jump at least " << setprecision(0)
         << STRONG_FACTOR << " x limit)" << endl;
    cout << "      weak                     : " << n_weak
         << "   (small jump and no flag: check by eye)" << endl;
    cout << "  Receiver loss-of-lock flags  : "
         << (n_flag_at_jump + n_flag_at_start + n_flag_only) << endl;
    cout << "      at an LG jump            : " << n_flag_at_jump << endl;
    cout << "      at the start of an arc   : " << n_flag_at_start
         << "   (normal: the receiver has just locked on)" << endl;
    cout << "      with no LG jump          : " << n_flag_only << endl;

    // ---- how to read it ----
    cout << endl << " HOW TO READ THIS" << endl
         << "  LG jump : change of LG from the previous epoch. Real slips are"
         << endl
         << "            19.0 cm (L1 only), 24.4 cm (L2 only), 5.4 cm (one"
         << " cycle on both)." << endl
         << "  limit   : the largest change accepted as normal for that time"
         << " step." << endl
         << "  x limit : how many times the limit the jump was." << endl
         << "  flag    : the receiver's own loss-of-lock flag, and on which"
         << " frequency." << endl
         << "  STRONG  : treat as a real slip." << endl
         << "  WEAK    : just over the limit and no flag. Often a false"
         << " alarm soon" << endl
         << "            after the satellite rises. Plot LG before trusting"
         << " it." << endl
         << line << endl;

    return 0;
}
