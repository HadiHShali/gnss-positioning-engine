// ============================================================================
//  rinex_obs.h
// ----------------------------------------------------------------------------
//  Data structures and function declarations for reading RINEX 3 observation
//  files (GPS satellites only).
//
//  A header answers "WHAT exists?". The matching rinex_obs.cpp answers
//  "HOW does it work?". Read the sections in order; each one builds on the
//  one before it.
//
//      Section 1   include guard + includes
//      Section 2   ObsRecord     (struct)       one satellite at one epoch
//      Section 3   ObsEpoch      (struct)       one epoch = time + satellites
//      Section 4   obsColumn     (inline func)  where an observable starts
//      Section 5   findObsIndex  (function)     which position holds "L1C"?
//      Section 6   parser entry points (functions)
//      Section 7   secondsBetween (inline func)  week-safe time difference
// ============================================================================


// ============================================================================
//  SECTION 1 -- INCLUDE GUARD AND INCLUDES
// ----------------------------------------------------------------------------
//  #pragma once : the compiler reads this file only once per .cpp, even if
//                 several headers include it. Without it the structs below
//                 would be defined twice -> "redefinition" error.
//  Includes     : only what THIS header uses itself:
//                 <string> for std::string, <vector> for std::vector.
// ============================================================================
#pragma once
#include <string>
#include <vector>


// ============================================================================
//  SECTION 2 -- ObsRecord : ONE SATELLITE AT ONE EPOCH
// ----------------------------------------------------------------------------
//  Why a struct : these are plain values that travel together. No rule ties
//                 one member to another, so everything is public.
//  Defaults     : every member has "= value", so a fresh ObsRecord means
//                 "nothing observed yet" instead of holding garbage.
//  Units        : pseudorange is in METERS, carrier phase is in CYCLES.
//                 Multiply phase by the wavelength (lambda = c / f) before
//                 combining it with anything in meters.
//  has_* flags  : RINEX writes a missing observation as blanks or 0.0, so
//                 the value alone cannot tell you it is missing. The flag can.
//  LLI          : bit 0 (value 1) = lock was lost, cycle slip possible
//                 bit 1 (value 2) = half-cycle ambiguity possible
//                 Test the slip bit on EACH frequency:
//                 (lli_l1 & 1) for L1 and (lli_l2 & 1) for L2.
// ============================================================================
struct ObsRecord {
    int    prn = 0;                 // GPS satellite number: 6 for "G06"

    // ---- pseudorange (code), in meters ----
    double pseudorange    = 0.0;    // C1W: P(Y)-code range on L1, meters
    double pseudorange_c2 = 0.0;    // C2W: P(Y)-code range on L2, meters
    bool   has_c2 = false;          // true if C2W was present this epoch

    // ---- carrier phase, in cycles (NOT meters) ----
    double phase_l1 = 0.0;          // L1C: carrier phase on L1, cycles
    double phase_l2 = 0.0;          // L2W: carrier phase on L2, cycles
    bool   has_phase_l1 = false;    // true if L1C was present this epoch
    bool   has_phase_l2 = false;    // true if L2W was present this epoch

    // ---- loss-of-lock indicator, one per phase ----
    int    lli_l1 = 0;              // LLI digit that follows L1C (0 = OK)
    int    lli_l2 = 0;              // LLI digit that follows L2W (0 = OK)
};


// ============================================================================
//  SECTION 3 -- ObsEpoch : ONE EPOCH = ONE TIME TAG + ALL ITS SATELLITES
// ----------------------------------------------------------------------------
//  Why a struct : same reason as ObsRecord -- plain data, no rule to protect.
//  Composition  : it holds a std::vector<ObsRecord>, one entry per GPS
//                 satellite seen at this epoch. A struct can contain other
//                 structs; that is how the file's shape (epoch -> satellites)
//                 becomes the data's shape.
//  Time, twice  : the calendar fields are kept as written in the file (for
//                 printing). gps_week + t_gps are for arithmetic (time
//                 differences, gap detection).
// ============================================================================
struct ObsEpoch {
    // ---- calendar time, exactly as written in the epoch line ----
    int    year = 0, month = 0, day = 0, hour = 0, minute = 0;
    double second = 0.0;            // seconds can have decimals

    // ---- the same instant in GPS time, for arithmetic ----
    double t_gps    = 0.0;          // seconds of GPS week (0 = Sunday 00:00:00)
    int    gps_week = 0;            // weeks since 1980-01-06

    int    epoch_flag = 0;          // 0 = OK, 1 = power failure since the
                                    // previous epoch, >1 = special event

    std::vector<ObsRecord> gps_records;   // one ObsRecord per GPS satellite
};


// ============================================================================
//  SECTION 4 -- obsColumn : WHERE DOES OBSERVABLE NUMBER i START IN A LINE?
// ----------------------------------------------------------------------------
//  Why a function : a pure calculation -- number in, number out, nothing to
//                   remember between calls.
//  Why "inline"   : the body lives in a header that many .cpp files include.
//                   "inline" lets each of them carry the same definition
//                   without a "multiple definition" linker error. Use this
//                   only for short, trivial functions.
//  The layout     : a RINEX 3 satellite line is
//                        3 characters   satellite ID ("G06")
//                     + 16 characters   per observable:
//                           14 value (F14.3) + 1 LLI + 1 signal strength
//                   For an observable that starts at column c:
//                        value = columns c .. c+13
//                        LLI   = column  c+14
//                        SSI   = column  c+15
//  Limits         : RINEX 3 only. "index" counts from 0.
// ============================================================================
inline int obsColumn(int index) {
    return 3 + 16 * index;          // skip the sat ID, then 16 chars per slot
}


// ============================================================================
//  SECTION 5 -- findObsIndex : WHICH POSITION HOLDS A GIVEN OBSERVABLE CODE?
// ----------------------------------------------------------------------------
//  Why a function   : a lookup with no memory of its own.
//  Declaration only : the line ends with ";" and has no body. The body goes
//                     in rinex_obs.cpp. Longer functions belong there, so
//                     that changing them does not recompile every file that
//                     includes this header.
//  Input            : obs_types = the codes read from the header record
//                     "SYS / # / OBS TYPES", in file order, for example
//                     ["C1C", "L1C", "D1C", "S1C", "C1W", "S1W", "C2W", ...]
//                     code      = the one you want, for example "L2W"
//  Output           : its position in the list (0-based), or -1 if absent.
//                     ALWAYS check for -1 before passing it to obsColumn().
//  const &          : both arguments are passed by const reference --
//                     no copy is made, and the function cannot modify them.
// ============================================================================
int findObsIndex(const std::vector<std::string>& obs_types,
                 const std::string& code);


// ============================================================================
//  SECTION 6 -- PARSER ENTRY POINTS : THE TWO WAYS TO READ A FILE
// ----------------------------------------------------------------------------
//  Why functions    : file name in, data out. Nothing is remembered after
//                     the call returns, so no class is needed.
//  Declaration only : the bodies are in rinex_obs.cpp.
//  Return by value  : returning a struct or a vector does not copy it element
//                     by element. The compiler builds it directly in the
//                     caller's variable, or moves it there.
//  These two functions are what the rest of the project calls. Everything
//  above exists to describe what they return.
// ============================================================================

// Read ONE epoch: the one at target_t_gps (seconds of GPS week),
// or the first epoch in the file if target_t_gps is 0.
ObsEpoch parseRinexObsOneEpoch(const std::string& filename,
                               double target_t_gps);

// Read EVERY epoch in the file, in file order.
std::vector<ObsEpoch> parseRinexAllEpochs(const std::string& filename);


// ============================================================================
//  SECTION 7 -- secondsBetween : A TIME DIFFERENCE THAT SURVIVES A WEEK CHANGE
// ----------------------------------------------------------------------------
//  Why it exists : t_gps restarts at 0 every Sunday. Subtracting two t_gps
//                  values across Saturday -> Sunday gives about -604800
//                  instead of the real step. Using gps_week as well fixes it.
//  Use it for    : every time step and every gap test. The same code then
//                  works for one day, for several days in a row, and for a
//                  run that crosses a week or a year.
//  Returns       : seconds from "earlier" to "later". The result is negative
//                  if the two epochs are passed in the other order.
//  Why inline    : one line, defined in the header (same reason as obsColumn).
// ============================================================================
inline double secondsBetween(const ObsEpoch& earlier, const ObsEpoch& later) {
    const double SECONDS_PER_WEEK = 604800.0;           // 7 * 86400
    return (later.gps_week - earlier.gps_week) * SECONDS_PER_WEEK   // weeks
         + (later.t_gps - earlier.t_gps);                // plus within the week
}
