//---------------------------------------------//
// Day 2: Haversine Distance Calculator
//---------------------------------------------//

//
// SECTION 1: Headers
#include <iostream>
#include <cmath>

using namespace std;

//
// SECTIPN 2: Constants
const double C_Light = 299792458.0; //speed of light (m/s)
const double R_EARTH = 6371000.0; //Earth Mean Radius

//
// SECTION 3: Pseudorange Function
// Formula: rho = c * Delta_t
double pseudorange(double travel_time_sec)
{
	return C_Light * travel_time_sec;

}

//
// SECTION 4: Haversine Functions
// Formula: a = sin^2(dphi/2) + cos(phi1)*cos(phi2)*sin^2(dlam/2)
//			c = 2*asin(sqrt(a))
//			d = R *c
const double M_PI = 3.14159265358979323846;

double haversine(double lat1_deg, double lon1_deg,
				 double lat2_deg, double lon2_deg)
{
	// SubStep1: convert to radians
	double phi1 = lat1_deg * M_PI / 180.0;
	double phi2 = lat2_deg * M_PI / 180.0;
	double lam1 = lon1_deg * M_PI / 180.0;
	double lam2 = lon2_deg * M_PI / 180.0;

	// SubStep2: Compute differences
	double d_phi = phi2 - phi1; // delta latitude
	double d_lam = lam2 - lam1; // delta longitude

	//SubStep3: Haversine Intermediate value
	double a = sin(d_phi / 2) * sin(d_phi / 2) +
		cos(phi1) * cos(phi2) * sin(d_lam / 2) * sin(d_lam / 2);

	//SubStep4: Central Angle
	double c = 2.0 * asin(sqrt(a));

	//SubStep5: Arc Distance
	return R_EARTH * c;

}

//
// SECTION 5: Main Program
int main()
{
	// Coordinates
	double rec_lat = 35.6892, rec_lon = 51.3890; // Receiver is in Tehran, Iran
	double sat_lat = 36.0, sat_lon = 52.0; // Satellite ground track 

	// Pseudorange
	double rho = pseudorange(0.0673); // travel time: 67.3 mili second
	cout << "=== GPS Distance Calculator ===" << endl << endl;
	cout << "Pseudorange : " << rho / 1000.0 << " km" << endl;

	// Haversine Surface Distance
	double d = haversine(rec_lat, rec_lon, sat_lat, sat_lon);
	cout << "Surface Distance: " << d / 1000.0 << " km" << endl;

	// Approximate Slant range
	double sat_alt = 20200000.0; // in meter
	double slant_range = sqrt(sat_alt * sat_alt + d * d);
	cout << "Slant Range: " << slant_range / 1000.0 << " km" << endl;

	// 
}