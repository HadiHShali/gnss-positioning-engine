// KalmanFilter.cpp : Actial implementation (the kitchen)
#include "KalmanFilter.h"  // pull in our own declaration

// Constructor: ClassName::ClassName means "this function belongs to KalmanFilter"
KalmanFilter::KalmanFilter(double x0, double P0, double q, double r)
{
    x_ = Eigen::MatrixXd(1, 1);
    x_(0, 0) = x0;
    P_ = Eigen::MatrixXd(1, 1);
    P_(0, 0) = P0;
    F_ = Eigen::MatrixXd(1, 1);
    F_(0, 0) = 1.0;
    H_ = Eigen::MatrixXd(1, 1);
    H_(0, 0) = 1.0;
    Q_ = Eigen::MatrixXd(1, 1);
    Q_(0, 0) = q;
    R_ = Eigen::MatrixXd(1, 1);
    R_(0, 0) = r;
    I_ = Eigen::MatrixXd::Identity(1, 1);

}

void KalmanFilter::predict()
{
    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::update(double measurement)
{
    Eigen::MatrixXd z(1, 1);
    z(0, 0) = measurement;

    Eigen::MatrixXd K = P_ * H_.transpose() * (H_ * P_ * H_.transpose() + R_).inverse();
    x_ = x_ + K * (z - H_ * x_);
    P_ = (I_ - K * H_) * P_;

}

double KalmanFilter::getState() const
{
    return x_(0, 0);
}

double KalmanFilter::getCovariance() const
{
    return P_(0, 0);
}