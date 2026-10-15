// real_data_position.cpp - Compute receiver position from real RINEX data
#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "gps_ls_solver.h"        // Week 6 solver
#include "utils.h"                  // ecefToLla
#include "satellites.h"            // Week 5 — parseRinexNav, computeSatPosECEF
#include "rinex_obs.h"             // NEW today
const double C_LIGHT = 299792458.0;
const double OMEGA_EARTH = 7.2921151467e-5;
using namespace std;

int main() {
    cout << "=== REAL GPS POSITIONING ===" << endl << endl;

    // STEP 1: Parse RINEX NAV (broadcast ephemerides)
    string nav_file = "../../data/BRDC00IGS_R_20240070000_01D_MN.rnx";
    vector<GpsEphemeris> ephs = parseRinexNav(nav_file);
    cout << "Parsed " << ephs.size() << " ephemerides" << endl;

    // STEP 2: Parse RINEX OBS (pseudoranges) Noon utc
    string obs_file = "../../data/BILL00USA_R_20240070000_01D_15S_MO.rnx";
    ObsEpoch epoch = parseRinexObsOneEpoch(obs_file, 43200.0);
    cout << "Epoch: " << epoch.year << "-" << epoch.month << "-" << epoch.day
        << " " << epoch.hour << ":" << epoch.minute << ":" << epoch.second << endl;
    cout << "GPS sec of week: " << epoch.t_gps << endl;
    cout << "Number of GPS observations: " << epoch.gps_records.size() << endl;

    if (epoch.gps_records.size() < 4) {
        cerr << "Need at least 4 satellites!" << endl;
        return 1;
    }


    // STEP 3: For each observed satellite, compute its position
    vector<Observation> obs;
    cout << endl << "Building observations..." << endl;
    // In real_data_position.cpp:

    for (const auto& rec : epoch.gps_records) {
        const GpsEphemeris* eph = findBestEphemeris(ephs, rec.prn, epoch.t_gps);
        if (!eph) continue;

        SatPosition sp = computeSatPosECEF(*eph, epoch.t_gps);

        // ── EARTH ROTATION CORRECTION (Sagnac) — ADD THIS ──
        double travel_time = rec.pseudorange / C_LIGHT;
        double rotation_angle = OMEGA_EARTH * travel_time;
        double cos_r = cos(rotation_angle);
        double sin_r = sin(rotation_angle);
        double sx_new = cos_r * sp.X + sin_r * sp.Y;
        double sy_new = -sin_r * sp.X + cos_r * sp.Y;
        sp.X = sx_new;
        sp.Y = sy_new;
        // ───────────────────────────────────────────────────

        double sat_clock_correction = eph->clk_bias
            + eph->clk_drift * (epoch.t_gps - eph->toe);

        Observation o;
        o.prn = rec.prn;
        o.sat_x = sp.X;
        o.sat_y = sp.Y;
        o.sat_z = sp.Z;
        o.pseudorange = rec.pseudorange + C_LIGHT * sat_clock_correction;
        obs.push_back(o);
    }
    cout << "Built " << obs.size() << " observation tuples" << endl << endl;

    if (obs.size() < 4) { cerr << "Not enough!" << endl; return 1; }

    // STEP 4: Solve!
    PositionSolution sol = solveGpsPosition(obs, false);
    if (!sol.converged) {
        cerr << "Did not converge!" << endl;
        return 1;
    }

    // STEP 5: Print results
    cout << "=== POSITION SOLUTION ===" << endl;
    cout << fixed << setprecision(3);
    cout << "Iterations: " << sol.iterations << endl;
    cout << "ECEF:" << endl;
    cout << "  X = " << sol.x << " m" << endl;
    cout << "  Y = " << sol.y << " m" << endl;
    cout << "  Z = " << sol.z << " m" << endl;

    double lat, lon, h;
    ecefToLla(sol.x, sol.y, sol.z, lat, lon, h);
    cout << endl << "LLA:" << endl;
    cout << setprecision(6);
    cout << "  Latitude:  " << lat << " deg" << endl;
    cout << "  Longitude: " << lon << " deg" << endl;
    cout << setprecision(2);
    cout << "  Height:    " << h << " m" << endl;

    cout << endl << "Clock bias: " << sol.clock_bias * 1e6 << " microseconds" << endl;
    cout << "RMS residual: " << sol.rms_residual << " m" << endl;
    cout << "Sigma-0:      " << sol.sigma0 << " m" << endl;

    return 0;
}