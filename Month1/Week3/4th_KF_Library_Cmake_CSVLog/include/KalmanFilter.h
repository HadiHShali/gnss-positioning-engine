// KalmanFilter .h : declares the interface (ther menu)

#pragma once // prevents file being included twice in the same build

#include <Eigen/Dense>

class KalmanFilter
{
public:
	// public interface: what users of this class allowed do
	// Everything written below this, can be accessed from outside the class.

	// Constroctor: a special function automatically called when an object is created: It initializes the filter.
	KalmanFilter(double x0, double P0, double q, double r);

	// Prediction
	void predict();

	// Update
	void update(double measurement);

	// Getter
	// getState(): This returns the current estimated state.
	double getState() const; // const: does not modify the object
	// getCovariance(): returns the current uncertainty estimate.
	double getCovariance() const;

private:
	// Everything below is hidden from (users) outside code.
	// This is encapsulation.
	// This protects the filter from accidental misuse.
	// The "_" suffix is a naming convention meaning: private member variable
	Eigen::MatrixXd x_;  // state estimate
	Eigen::MatrixXd P_;   // covariance
	Eigen::MatrixXd F_;   // state transition
	Eigen::MatrixXd H_;   // measurement model
	Eigen::MatrixXd Q_;   // process noise
	Eigen::MatrixXd R_;   // measurement noise
	Eigen::MatrixXd I_;   // identity matrix


};

