#pragma once
#include <vector>
#include <string>

struct GpsEphemeris {
    int prn;
    int year, month, day, hour, minute;
    double second;
    double clk_bias, clk_drift, clk_drift_rate;
    double IODE, Crs, delta_n, M0;
    double Cuc, e, Cus, sqrt_a;
    double toe, Cic, Omega0, Cis;
    double i0, Crc, omega, Omega_dot;
    double i_dot, L2_codes, gps_week, L2P_flag;
    double sv_accuracy, sv_health, TGD, IODC;
    double trans_time, fit_interval;
};

struct SatPosition { double X, Y, Z; };

// Public API
std::vector<GpsEphemeris> parseRinexNav(const std::string& filename);

SatPosition computeSatPosECEF(const GpsEphemeris& eph, double t_gps);

const GpsEphemeris* findBestEphemeris(
    const std::vector<GpsEphemeris>& ephs, int prn, double t_gps);
	
// --------------------------------------------------------------------------------------------------------

double computeElevation(double rec_x, double rec_y, double rec_z,
                        double sat_x, double sat_y, double sat_z);
						
double computeAzimuth (double rec_x, double rec_y, double rec_z, double sat_x, double sat_y, double sat_z);

// ----------------------------------------------------------------------------------------------------------

// Klobuchar Ionospheric Correction parameters
// Read from Rinex Nav header (GPSA + GPSB lines)
struct KlobucharParams
{
    double alpha[4] = {};   // amplitude coefficients a0..a3 [s, s/sc, s/sc^2, s/sc^3]  (sc = semicircle)
    double beta[4] = {};   // period coefficients    b0..b3 [s, s/sc, s/sc^2, s/sc^3]
    bool valid = false;     // true only if BOTH alpha and beta were found
};

// Function to parse Klobuchar params from the RINEX Nav file header
KlobucharParams parseKlobucharFromNav(const std::string& filename);

// compute Ionospheric delay in meters using Klobuchar model
// input: KlobucharParams, receiver latitude, longitude, satellite azimuth, elevation, and GPS seconds of week
double klobucharIonoDelay (const KlobucharParams& kp, double rec_lat_deg, double rec_lon_deg,double sat_az_deg,double sat_el_deg, double t_gps_sec);

// --------------------------------------------------------------------------------------------------------------------------------------------
// Saastamoinen Tropospheric Delay (returns delays in meters)
// it depends only to receiver latitude, height and satellite elevation angle.
double saastamoinenTropoDelay(double rec_lat_deg, double rec_height_m, double sat_el_deg);
