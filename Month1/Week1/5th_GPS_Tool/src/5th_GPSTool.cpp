//==============================//
//   DAY 5 : Complete GPS Tool  //
//==============================//

// SECTION 1 : Headers
#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
using namespace std;

// SECTION 2: Constant
const double WGS84_A = 6378137.0;
const double WGS84_E2 = 0.00669437999014;
const double C_Light = 299792458.0;
const double Mask_Angl = 10.0;

// SECTION 3: Unit Conversion
double toRadians(double deg)
{
	return deg * M_PI / 180.0;
}
double toDegrees(double rad)
{
	return rad * 180.0 / M_PI;
}

//SECTION 4: Pseudorange
double pseudorange(double t)
{
	return t * C_Light;
}

// SECTION 5: LLA to ECEF
double Nv(double phi)
{
	return WGS84_A / sqrt(1.0 - WGS84_E2 * sin(phi) * sin(phi));
}

void lla2ecef(double lat, double lon, double alt, double& X, double& Y, double& Z)
{
	double phi = toRadians(lat), lam = toRadians(lon), N = Nv(phi);
	X = (N + alt) * cos(phi) * cos(lam);
	Y = (N + alt) * cos(phi) * sin(lam);
	Z = (N * (1.0 - WGS84_E2) + alt) * sin(phi);
}

// SECTION 6: Sat Classes
class satellite
{
public:
	int prn;
	double X, Y, Z, rho, elevation, azimuth;
	bool visible;

	satellite(int p, double x, double y, double z, double pr)
		: prn(p), X(x), Y(y), Z(z), rho(pr), elevation(0), azimuth(0), visible(false) { };

	void computeElevAZ(double Xr, double Yr, double Zr, double lat, double lon)
	{
		double phi = toRadians(lat), lam = toRadians(lon);
		double dX = X - Xr, dY = Y - Yr, dZ = Z - Zr;
		double e = -sin(lam) * dX + cos(lam) * dY;
		double n = -sin(phi) * cos(lam) * dX - sin(phi) * sin(lam) * dY + cos(phi) * dZ;
		double u = cos(phi) * cos(lam) * dX + cos(phi) * sin(lam) * dY + sin(phi) * dZ;
		elevation = atan2(u, sqrt(e * e + n * n)) * 180.0 / M_PI;
		azimuth = atan2(e, n) * 180.0 / M_PI;
		if (azimuth < 0) azimuth += 360;
		visible = (elevation >= Mask_Angl);
	}
	void print()
	{
		cout << "PRN: " << prn << " ,El: " << elevation << " ,AZ: " << azimuth
			<< " " << (visible ?"[VISIBLE]":"[MASKED]") << endl;
	}
};

// SECTION 7: MAin Program
int main()
{
	cout << "===GPS Receiver Simulation===" << endl << endl;
	double rec_lat = 35.6892; 
	double rec_lon = 51.3890;
	double rec_alt = 1191.0;
	double Xr, Yr, Zr;
	lla2ecef(rec_lat, rec_lon, rec_alt, Xr, Yr, Zr);
	cout << "Receiver ECEF: Xr= " << Xr << ",Yr= " << Yr << ",Zr= " << Zr << endl << endl;

	satellite sats[4] = {
		satellite(1,15600000, 7540000,20140000,pseudorange(0.0673)),
		satellite(2,-9000000,20000000,15000000,pseudorange(0.0703)),
		satellite(3,21000000, 3500000,14000000,pseudorange(0.0743)),
		satellite(4, 5000000,-19000000,18000000,pseudorange(0.0167)),
	};
	int vis = 0;
	for (int i = 0; i < 4; i++)
	{
		sats[i].computeElevAZ(Xr, Yr, Zr, rec_lat, rec_lon);
		sats[i].print();
		if (sats[i].visible) vis++;
	}
	cout << endl << "Visible: " << vis << "/4 Fix: " << (vis > 4 ? "POSSIBLE" : "NOT POSSIBLE") 
		<< endl;
	return 0;
}

