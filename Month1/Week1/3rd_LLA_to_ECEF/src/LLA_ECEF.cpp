//
// Day 3: LLA to ECEF converter
//

//SECTION 1: Headers
#define _USE_MATH_DEFINES // required on windows visual studio for M_PI which is pi=3.1415
#include <iostream>
#include <cmath>

using namespace std; // instead of wrting "std::cout" we can write "cout"

//SECTION 2: WGS84 constants
const double WGS84_A = 6378137.0;
const double WGS84_E2 = 0.00669437999014;

// SECTION 3: Prime Vertical Radius = N, it is a function of phi
// Formula: N = a/sqrt(1-e2*sin^2(phi))
double PrimeVerticalRadius(double phi_rad)
{
	return WGS84_A / sqrt(1 - WGS84_E2 * sin(phi_rad) * sin(phi_rad));
}

//SECTION 4: LLA to ECEF
// Formula:
// X = (N + h) * cos(phi) * cos(lambda)
// Y = (N + h) * cos(phi) * sin(lambda)
// Z = (N * (1-e2) + h) * sin(phi)
void lla2ecef(double lat_deg, double lon_deg, double alt_m, double &X, double &Y, double &Z) // &X, &Y, and &Z are the output 
{
	double phi = lat_deg * M_PI / 180.0;
	double lam = lon_deg * M_PI / 180.0;
	double N = PrimeVerticalRadius(phi);
	X = (N + alt_m) * cos(phi) * cos(lam);
	Y = (N + alt_m) * cos(phi) * sin(lam);
	Z = (N * (1.0 - WGS84_E2) + alt_m) * sin(phi);
}


//SECTION 5: ECEF to LLA (Bowring iteration)
// Longitude is exact [lon = atan2(Y,X)]!!
// Latitude needs iterations to get converged, since it depends on the phi (phi appears on both sides)
// phi = atan2(Z + e2*N*sin(phi), p)

void ecef2lla(double X, double Y, double Z, double &lat_deg, double &lon_deg, double &alt_m)
{
	double p = sqrt(X * X + Y * Y);
	double lam = atan2(Y, X);
	double phi = atan2(Z, p * (1.0 - WGS84_E2)); // Initial estimate
	for (int i = 0; i < 6; i++)
	{
		double N = PrimeVerticalRadius(phi);
		phi = atan2(Z + WGS84_E2 * N * sin(phi), p);
	}
	double N = PrimeVerticalRadius(phi);
	lat_deg = phi * 180.0 / M_PI;
	lon_deg = lam * 180.0 / M_PI;
	alt_m = p / cos(phi) - N;

}

// SECTION 6: Main program
int main()
{
	// lat lon and altitude:
	double lat = 35.6892, lon = 51.3890, alt = 1191.0;
	// X, Y , Z
	double X, Y, Z;

	//Forward: LLA > ECEF
	lla2ecef(lat, lon, alt, X, Y, Z);
	cout << "Input (LLA): lat=" << lat << " lon=" << lon << " alt=" << alt << endl;
	cout << "Output (ECEF): X=" << X << " Y=" << Y << " Z=" << Z << endl << endl;

	//Inverse: ECEF > LLA (round trip check)
	double lat2, lon2, alt2;
	ecef2lla(X, Y, Z, lat2, lon2, alt2);
	cout << "ROUND-TRIP: lat=" << lat2 << " lon=" << lon2 << " alt=" << alt2 << endl << endl;

	// Additional: 4 satellites positions (ECEF > LLA)
	double sats[4][3] = { {15600000,7540000,20140000},{-9000000,20000000,15000000},{21000000,3500000,14000000},
		{5000000,-19000000,18000000} };

	for (int i = 0; i < 4; i++)
	{
		double sat_lat, sat_lon, sat_alt;
		ecef2lla(sats[i][0], sats[i][1], sats[i][2], sat_lat, sat_lon, sat_alt);
		cout << "sat " << i + 1 << ": lat=" << sat_lat << " lon=" << sat_lon << " alt=" << sat_alt / 1000.0 << "km" << endl;
	}
	return 0;
}