// app_kf_1d.cpp : uses KalmanFilter library

#include <iostream>
#include "KalmanFilter.h" //pull in our library
using namespace std;

int main()
{
	// create a Kalman Filter same as Werek 2 Day 1
	KalmanFilter kf(1150.0, 500.0, 1.0, 225.0);

	double gps[] = { 1205.3, 1178.9, 1198.2, 1183.7, 1204.1,
					1176.5, 1195.8, 1188.3, 1201.9, 1185.4 };

	cout << "=== Modular KF Library Test===" << endl;
	cout << "Epoch  GPS(m)    KF_Est(m)   Uncertainty" << endl;

	for (int i = 0; i < 10; i++)
	{
		kf.predict();
		kf.update(gps[i]);
		cout << "  " << i + 1
			<< "    " << gps[i]
			<< "    " << kf.getState()
			<< "    " << kf.getCovariance() << endl;
	}

	return 0;

}