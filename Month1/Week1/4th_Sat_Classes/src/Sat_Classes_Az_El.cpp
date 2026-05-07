//
// Day 4: GPS Satellites Class, Azimuth and Elevation Angles
//

// SECTION 1: Headers
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
using namespace std;

// SECTION 2: Constants
const double WGS84_A = 6378137.0;
const double WGS84_E2 = 0.00669437999014;
const double Mask_Deg = 10.0;

// SECTION 3: Satellite Classes
class Satellite
{
public:

	int prn;
	double X, Y, Z;
	double pseudorange;
	double elevation;
	double azimuth;
	bool visible;

	// Constructor: the name of the function should be the same as the name of the calss
	// that's how C++ knows it's a constructor.
	Satellite(int prn_, double x, double y, double z, double pr) //These are the inputs you pass when creating the object.
		: prn(prn_), X(x), Y(y), Z(z), pseudorange(pr), elevation(0), azimuth(0), visible(false) 
	{}
	// This is the initializer list — the : starts it. It directly initializes each member variable:
	// for example prn(prn_): set prn = the value passed as prn_
	// "{}" at the end is the body of the constructor. it is necessary. We can NOT remove it

	// SubSsection 1: Compute Elevation anf Azimuth Angles
	void computeElevAz(double Xr, double Yr, double Zr, double rec_lat_deg, double rec_lon_deg) 
	{
		double phi = rec_lat_deg * M_PI / 180.0;
		double lam = rec_lon_deg * M_PI / 180.0;
		double dX = X - Xr, dY = Y - Yr, dZ = Z - Zr;
		double e = -sin(lam) * dX + cos(lam) * dY;
		double n = -sin(phi) * cos(lam) * dX - sin(phi) * sin(lam) * dY + cos(phi) * dZ;
		double u = cos(phi) * cos(lam) * dX + cos(phi) * sin(lam) * dY + sin(phi) * dZ;
		elevation = atan2(u, sqrt(e * e + n * n)) * 180.0 / M_PI;
		azimuth = atan2(e, n) * 180.0 / M_PI;
		if (azimuth < 0) azimuth += 360.0;
		visible = (elevation >= Mask_Deg);
	}

	// SubSection 2: Print Report
	void print() {
		cout << "PRN " << prn
			<< "  El=" << elevation << "deg"
			<< "  Az=" << azimuth << "deg"
			<< "  Range=" << pseudorange / 1000.0 << " km"
			<< (visible ? "  [VISIBLE]" : "  [MASKED]") << endl;
	}

};

// SECTION 4: MAin program
int main() {
	double Xr = 3352834.0, Yr = 4076956.0, Zr = 3691542.0;
	double rec_lat = 35.6892, rec_lon = 51.3890;

	// SubSection1 — CREATE SATELLITE OBJECTS
	Satellite sats[4] = {
		Satellite(1, 15600000,  7540000, 20140000, 20200000),
		Satellite(2, -9000000, 20000000, 15000000, 21100000),
		Satellite(3, 21000000,  3500000, 14000000, 22300000),
		Satellite(4,  5000000,-19000000, 18000000,  5000000),
	};

	// SubSection2 — VISIBILITY REPORT
	cout << "=== GPS Satellite Visibility Report ===" << endl;
	int vis = 0;
	for (int i = 0; i < 4; i++) {
		sats[i].computeElevAz(Xr, Yr, Zr, rec_lat, rec_lon);
		sats[i].print();
		if (sats[i].visible) vis++;
	}
	cout << endl << "Visible: " << vis << "/4  Fix possible: "
		<< (vis >= 4 ? "YES" : "NO") << endl;
	return 0;
}
