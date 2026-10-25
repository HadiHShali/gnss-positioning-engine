#pragma once
// interpolation_crctd.h - corrected copy of interpolation.h (Month 3 / Week 14 review, 2026-10-24)
// The original header stays untouched; only spp_processor_precise_crctd.cpp includes this one.
// Changes are marked [CRCTD]:
//   1. interpolateSatPosition() EXTRAPOLATED past the last SP3 epoch. A daily IGS SP3 file
//      ends at 23:45, so every epoch from 23:45 to 23:59:45 was extrapolated with a
//      9th-degree polynomial: up to 1.6 m of orbit error after 15 min (synthetic test),
//      vs 0.5 mm in the middle of the day. Now it returns false outside the SP3 time span
//      and the processor falls back to broadcast orbit AND clock there.
//   2. SP3 marks a missing position as 0.000000 0.000000 0.000000. Those records were
//      used as real nodes (a satellite "at the Earth's centre"). Now they are skipped.
//   3. Both functions interpolated ACROSS data gaps without limit (e.g. a satellite missing
//      for 2 h). Now: orbit -> false if the two SP3 nodes around t are more than
//      MAX_SP3_GAP_S apart; clock -> false if the two CLK epochs are more than
//      MAX_CLK_GAP_S apart. A clock is not a straight line over hours (ns = metres).
// Unchanged: lagrangeInterp(), linearClockInterp(), the 10-node window choice.
#include <vector>
#include <string>
#include <cmath>
#include "sp3_reader.h"
#include "clk_reader.h"

const double MAX_SP3_GAP_S = 1200.0;   // [CRCTD] 15-min SP3: one node apart = 900 s; allow some slack
const double MAX_CLK_GAP_S = 300.0;    // [CRCTD] 30-s or 5-min CLK: anything longer is a gap

// 10th-order Lagrange interpolation for one coordinate
inline double lagrangeInterp(const std::vector<double>& times,
    const std::vector<double>& values,
    double t)
{
    double result = 0.0;
    int n = static_cast<int>(times.size());
    for (int i = 0; i < n; i++) {
        double term = values[i];
        for (int j = 0; j < n; j++) {
            if (j != i) term *= (t - times[j]) / (times[i] - times[j]);
        }
        result += term;
    }
    return result;
}

// Linear interpolation for clocks -- deliberately NOT Lagrange
inline double linearClockInterp(double t1, double clk1,
    double t2, double clk2, double t)
{
    return clk1 + (clk2 - clk1) * (t - t1) / (t2 - t1);
}


// Given a satellite's PRN and a target time, find that satellite's 10 nearest SP3 epochs,
// then use Lagrange interpolation to compute its ECEF position at target_t -> out_pos[3].
// Returns false if the satellite has < 10 valid nodes, if target_t is outside the
// SP3 time span [CRCTD], or if target_t falls in a data gap [CRCTD].
inline bool interpolateSatPosition(const std::vector<Sp3Epoch>& epochs,
    const std::string& prn,
    double target_t,
    double out_pos[3])
{

    // Phase 1 — Collect This Satellite's Data Across All Epochs
    std::vector<double> times;
    std::vector<double> xs, ys, zs;
    for (const auto& ep : epochs) {
        auto it = ep.sat_pos.find(prn);
        if (it != ep.sat_pos.end()) {
            // [CRCTD] skip SP3 "missing" records (all three coordinates 0.000000)
            if (it->second[0] == 0.0 && it->second[1] == 0.0 && it->second[2] == 0.0) continue;
            times.push_back(ep.t_gps);
            xs.push_back(it->second[0]);
            ys.push_back(it->second[1]);
            zs.push_back(it->second[2]);
        }
    }
    if (times.size() < 10) return false;   // not enough data for order-10

    // [CRCTD] interpolate only, never extrapolate: t must lie between the first and last node
    if (target_t < times.front() || target_t > times.back()) return false;

    // [CRCTD] the two nodes around t must not be separated by a data gap
    for (size_t i = 0; i + 1 < times.size(); i++) {
        if (times[i] <= target_t && target_t <= times[i + 1]) {
            if (times[i + 1] - times[i] > MAX_SP3_GAP_S) return false;
            break;
        }
    }

    // Phase 2: Find the 10 epochs closest to target_t (5 before, 5 after ideally)
    // Simplified: find nearest index, take a centered window of 10
    int nearest = 0;
    double best_diff = 1e18;
    for (size_t i = 0; i < times.size(); i++) {
        double diff = std::abs(times[i] - target_t);
        if (diff < best_diff) { best_diff = diff; nearest = static_cast<int>(i); }
    }
    int start = nearest - 5;
    if (start < 0) start = 0;
    if (start + 10 > static_cast<int>(times.size()))
        start = static_cast<int>(times.size()) - 10;

    std::vector<double> wt(times.begin() + start, times.begin() + start + 10);
    std::vector<double> wx(xs.begin() + start, xs.begin() + start + 10);
    std::vector<double> wy(ys.begin() + start, ys.begin() + start + 10);
    std::vector<double> wz(zs.begin() + start, zs.begin() + start + 10);

    out_pos[0] = lagrangeInterp(wt, wx, target_t);
    out_pos[1] = lagrangeInterp(wt, wy, target_t);
    out_pos[2] = lagrangeInterp(wt, wz, target_t);
    return true;
}



// Linear-interpolate a satellite's clock bias at target_t.
// Returns false if the satellite isn't found, target_t is outside the available
// CLK epoch range, or [CRCTD] the two bracketing epochs are a data gap apart.
inline bool interpolateSatClock(const std::vector<ClkEpoch>& epochs,
                                  const std::string& prn,
                                  double target_t,
                                  double& out_clock)
{
    // Find the two epochs bracketing target_t for this satellite
    double t1 = -1, t2 = -1, c1 = 0, c2 = 0;
    for (const auto& ep : epochs) {
        auto it = ep.sat_clock.find(prn);
        if (it == ep.sat_clock.end()) continue;
        if (ep.t_gps <= target_t) { t1 = ep.t_gps; c1 = it->second; }
        if (ep.t_gps >= target_t && t2 < 0) { t2 = ep.t_gps; c2 = it->second; }
    }
    if (t1 < 0 || t2 < 0) return false;
    if (t1 == t2) { out_clock = c1; return true; }   // exact match
    if (t2 - t1 > MAX_CLK_GAP_S) return false;      // [CRCTD] no straight line across a gap
    out_clock = linearClockInterp(t1, c1, t2, c2, target_t);
    return true;
}
