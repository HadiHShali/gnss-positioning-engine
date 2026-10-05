#pragma once
#include <vector>
#include "rinex_obs.h"

constexpr double C_LIGHT = 299792458.0;
constexpr double GPS_F1 = 1575.42e6;
constexpr double GPS_F2 = 1227.60e6;
constexpr double LAMBDA_L1 = C_LIGHT / GPS_F1;   // ~0.190294 m
constexpr double LAMBDA_L2 = C_LIGHT / GPS_F2;   // ~0.244210 m

// Geometry-free combination: LG = phase_l1 - phase_l2 (in CYCLES).
// Per Day 1: smooth ionosphere variation normally; a SUDDEN jump signals a likely cycle slip.

// Geometry-free combination, in METERS (not raw cycles!)
inline double geometryFree(double phase_l1_cycles, double phase_l2_cycles) {
    return LAMBDA_L1 * phase_l1_cycles - LAMBDA_L2 * phase_l2_cycles;
}

// Simple slip flag: compare consecutive-epoch LG values for one satellite.
// threshold_cycles: start with ~1.0 cycle as a conservative first pass
// (Day 1 noted single L1 slips are only ~0.221 LG cycles -- expect to
// tune this down once you see real data behavior).

//Why this works as a slip detector : per Day 1's theory, the ionosphere (which is what LG isolates)
// changes smoothly over time — normal epoch-to-epoch drift is tiny. A sudden large jump breaks that
// smooth pattern, which is exactly what happens when N (the integer ambiguity) abruptly changes due to a lost lock.

//Why 1.0 is called "conservative" in the comment : the docstring itself flags a tension worth noticing
// — Day 1's research found single L1-only slips show up as only ~0.221 LG cycles, much smaller than this 1.0 threshold.
// That means, as written, this function would miss genuine single-cycle L1 slips (0.221 < 1.0, 
// so likelySlip returns false even when a real slip occurred). It's deliberately set loose for a first pass — 
// catching only larger, more obvious slips — with the expectation you'll lower it once you see how noisy real LG data
// actually looks and can tell true slips apart from normal noise.

inline bool likelySlip(double lg_prev, double lg_curr, double threshold_cycles = 1.0)
{
	return std::abs(lg_curr - lg_prev) > threshold_cycles;
}


// Melbourne-Wubbena combination -- second, independent slip detector
// Needs C1W/C2W pseudoranges too (not just phase) -- p1_m and p2_m are your existing rec.pseudorange / rec.pseudorange_c2 fields from Week 13.
inline double melbourneWubbena(
    double p1_m,
    double p2_m,
    double phase_l1_cycles,
    double phase_l2_cycles)
{
    const double lambda_WL =
        C_LIGHT / (GPS_F1 - GPS_F2);

    const double L_WL =
        lambda_WL *
        (phase_l1_cycles - phase_l2_cycles);

    const double P_NL =
        (GPS_F1 * p1_m + GPS_F2 * p2_m)
        / (GPS_F1 + GPS_F2);

    return L_WL - P_NL;
}
