// ============================================================================
//  cycle_slip.h
// ----------------------------------------------------------------------------
//  Cycle-slip detection with the geometry-free (LG) phase combination.
//
//  Everything here is a constant or a short pure function, so it all lives
//  in the header and there is no cycle_slip.cpp.
//
//      Section 1   include guard + includes
//      Section 2   GPS signal constants   (constexpr)
//      Section 3   geometryFree           (inline function)  LG in meters
//      Section 4   likelySlip             (inline function)  did LG jump?
//      Section 5   slipThreshold          (inline function)  limit per time step
// ============================================================================


// ============================================================================
//  SECTION 1 -- INCLUDE GUARD AND INCLUDES
// ----------------------------------------------------------------------------
//  #pragma once : same job as in rinex_obs.h -- read this file once per .cpp.
//  <cmath>      : std::abs (used in Section 4) and std::exp (Section 5).
//  Not included : <vector> and "rinex_obs.h". Nothing in this header uses a
//                 vector, an ObsRecord or an ObsEpoch. The functions take
//                 plain doubles, so this header stands on its own.
// ============================================================================
#pragma once
#include <cmath>


// ============================================================================
//  SECTION 2 -- GPS SIGNAL CONSTANTS
// ----------------------------------------------------------------------------
//  Why constexpr  : these are fixed values known at compile time. constexpr
//                   makes them read-only and lets the compiler do the two
//                   divisions below once, while compiling.
//  Safe in a header : a constexpr variable is private to each .cpp that
//                   includes this file, so there is no "multiple definition"
//                   linker error and no "inline" is needed.
//  Wavelength     : lambda = c / f   (NOT 1 / f).
//  Values         : c is the exact SI value. f1 and f2 are the GPS L1 and L2
//                   carrier frequencies (IS-GPS-200).
// ============================================================================
constexpr double C_LIGHT   = 299792458.0;        // speed of light, m/s (exact)
constexpr double GPS_F1    = 1575.42e6;          // L1 carrier frequency, Hz
constexpr double GPS_F2    = 1227.60e6;          // L2 carrier frequency, Hz
constexpr double LAMBDA_L1 = C_LIGHT / GPS_F1;   // L1 wavelength, 0.190294 m
constexpr double LAMBDA_L2 = C_LIGHT / GPS_F2;   // L2 wavelength, 0.244210 m


// ============================================================================
//  SECTION 3 -- geometryFree : THE LG COMBINATION, IN METERS
// ----------------------------------------------------------------------------
//  Formula  : LG = lambda1 * phi1 - lambda2 * phi2
//  Input    : the two carrier phases in CYCLES, exactly as stored in
//             ObsRecord (phase_l1, phase_l2).
//  Output   : LG in METERS.
//  Why multiply first : phi1 - phi2 in raw cycles does NOT cancel the range.
//             It leaves 1.16 cycles per meter of range. Converting each
//             phase to meters first makes range, clocks and troposphere
//             cancel exactly.
//  What is left : ionosphere (smooth in time) + ambiguities (constant until
//             a slip). So a sudden jump in LG means a cycle slip.
//  Caller must check : has_phase_l1 AND has_phase_l2 are both true.
// ============================================================================
inline double geometryFree(double phase_l1_cycles, double phase_l2_cycles) {
    return LAMBDA_L1 * phase_l1_cycles      // L1 phase: cycles -> meters
         - LAMBDA_L2 * phase_l2_cycles;     // L2 phase: cycles -> meters
}


// ============================================================================
//  SECTION 4 -- likelySlip : DID LG JUMP BETWEEN TWO CONSECUTIVE EPOCHS?
// ----------------------------------------------------------------------------
//  Test      : |LG(now) - LG(previous)| > threshold
//  Units     : everything in METERS, because LG is in meters.
//  Threshold : default 0.05 m (5 cm). Slip sizes in LG:
//                  L1 only, 1 cycle          19.0 cm   detected
//                  L2 only, 1 cycle          24.4 cm   detected
//                  1 cycle on both            5.4 cm   barely detected
//              Navipedia's formula gives 4.9 cm for 15 s data, which is the
//              BILL00USA sampling interval.
//  Limits    : (1) Use it only for two epochs of the SAME satellite in the
//                  SAME arc. After a long gap, start a new arc instead of
//                  calling this.
//              (2) LG cannot see slip pairs near the ratio 77:60, for
//                  example 9 cycles on L1 with 7 on L2. The Melbourne-
//                  Wubbena detector covers those.
// ============================================================================
inline bool likelySlip(double lg_prev, double lg_curr,
                       double threshold_m = 0.05) {
    return std::abs(lg_curr - lg_prev) > threshold_m;   // all three in meters
}


// ============================================================================
//  SECTION 5 -- slipThreshold : HOW BIG A JUMP IS TOO BIG FOR THIS TIME STEP?
// ----------------------------------------------------------------------------
//  Why it exists : the ionosphere has more time to change over a longer time
//                  step, so the allowed change in LG must grow with the step.
//                  One fixed number cannot suit 1 s data, 30 s data and a
//                  comparison across a short dropout at the same time.
//  Formula       : threshold(dt) = a0 - (a0 / 2) * exp(-dt / T0)
//                      a0 = 1.5 * (lambda2 - lambda1) = 8.1 cm
//                      T0 = 60 s
//                  (ESA Navipedia, geometry-free cycle-slip detector)
//  Values        : dt =  1 s -> 4.1 cm        dt = 30 s -> 5.6 cm
//                  dt = 15 s -> 4.9 cm        dt = 60 s -> 6.6 cm
//  Input         : dt_seconds = time between the two epochs being compared.
//                  Get it from secondsBetween() in rinex_obs.h.
//  Output        : the threshold in METERS, ready to pass to likelySlip().
// ============================================================================
inline double slipThreshold(double dt_seconds) {
    const double A0 = 1.5 * (LAMBDA_L2 - LAMBDA_L1);    // 0.0809 m
    const double T0 = 60.0;                             // seconds
    return A0 - 0.5 * A0 * std::exp(-dt_seconds / T0);  // meters
}
