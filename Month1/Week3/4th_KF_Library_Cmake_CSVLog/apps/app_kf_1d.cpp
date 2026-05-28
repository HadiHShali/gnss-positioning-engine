// app_kf_1d.cpp : uses KalmanFilter library

#include <iostream>
#include <fstream>
#include "KalmanFilter.h" //pull in our library
using namespace std;

int main()
{
	// create a Kalman Filter same as Werek 2 Day 1
	KalmanFilter kf(1150.0, 500.0, 1.0, 225.0);

	double gps[] = { 1205.3, 1178.9, 1198.2, 1183.7, 1204.1,
					1176.5, 1195.8, 1188.3, 1201.9, 1185.4 };
	// Open the CSV file
	ofstream csvv("../../data/kf_1d_output.csv");  //csvv: is the variable name
	// got to data folder and save the file there
	
	// write the header now
	csvv << "epoch,gps_alt,kf_alt,uncertainty" << endl;

	// print to the screen too
	cout << "=== Modular KF Library with CSV Logging===" << endl;
	cout << "Epoch  GPS(m)    KF_Est(m)   Uncertainty" << endl;

	for (int i = 0; i < 10; i++)
	{
		kf.predict();
		kf.update(gps[i]);

		double est = kf.getState();
		double unc = kf.getCovariance();

		// Write one CSV row: 
		csvv << (i + 1) << ","
			<< gps[i] << ","
			<< est << ","
			<< unc << endl;

		// Also print on the screen
		cout << " " << (i + 1)
			<< "  " << gps[i]
			<< "  " << est
			<< "  " << unc << endl;
	}

	// Close the file
	csvv.close();

	cout << endl << "Results saved to kf_1d_output.csv" << endl;

	return 0;

}