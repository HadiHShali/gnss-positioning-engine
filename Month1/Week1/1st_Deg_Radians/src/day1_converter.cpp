//---------------------------------------------------------------//	
// Day 1: GPS coordinates converter (deg to rad and rad to deg).//
//--------------------------------------------------------------//
// 
// SECTION 1 : Headers
// 
// iostream : gives us "cout" and "cin" to print to the console and read from the console
// cmath : gives us "M_PI", "cos", "sin", "sqrt", "atan2", "pow" and other math functions - All GPS math uses
#include <iostream>
#include <cmath>
using namespace std; // this line allows us to write "cout" instead on "std::cout" and "cin" instead of "std::cin"

//
// SECTION 2: WGS84 constants
// 
// WGS84 = World Geodetic System 1984 - the mathematical model of the earth used by All GPS calculations.
// a = wquatorial radius (semi-major axis) in meters
// b = polar radius (semi-minor axis) in meters
// e2 = eccentricity squared - how much Earth bulges at the equator compared to the poles. Used in GPS math.
const double WGS84_A = 6378137.0; //Earth's equatorial radius 
const double WGS84_B = 6356752.314; //Earth's polar radius
const double WGS84_E2 = 1 - (WGS84_B * WGS84_B) / (WGS84_A * WGS84_A); //Eccentricity squared: 1-(b^2/a^2)=0.00669437999014

//
// SECTION 3: Core Conversion Functions
// 
// Formula: radians = degrees * (pi/180) >> This is called before every sin(), cos(), atan2(), etc calls in all GPS positioning code
// The output of the function is "double". Name of the funtion is "toRadians". input is a type "double" and called "deg" 
double M_PI = 3.14159265358979323846; // Define M_PI if not defined in cmath
double toRadians(double deg) { 
	return deg * (M_PI / 180.0); // 180.0 NOT 180 >> forces floating-point division
}
// Formula: degrees = radians * (180 / pi)
// This one is used when converting ECEF position output back to human-readable lat/lon coordinates
double toDegress(double rad) {
	return rad * (180.0 / M_PI);
}

// 
// SECTION 4: Main Program
//
int main() {
	// Real GPS coordinates: Tehran, Iran
	// same format as an NMEA GGA sentence from any GPS receiver. Example: $GPGGA,123519,35.6892,N,51.3890,E,1,08,0.9,545.4,M,46.9,M,,*47
	double lat = 35.6892; // degrees North
	double lon = 51.3890; // degrees East
	double alt = 1191.0; // meters above WGS84 elliposid (NOT sea level) 

	// Convert and Display the results
	double lat_radians = toRadians(lat);
	double lon_radians = toRadians(lon);
	cout << "=== GPS Coordinates Converter ===" << endl;
	cout << "Input Latitude: " << lat << " deg = " << lat_radians << " radians" << endl;
	cout << "Input Longitude: " << lon_radians << " deg = " << lon_radians << " radians" << endl;
	cout << "Altitude: " << alt << " meters" << endl;
	return 0; // Why do we need to return 0? In C++, the main function must return an integer. 
	// Returning 0 typically indicates that the program executed successfully without errors. 
	// It's a convention that signals to the operating system that the program finished as expected.

}

