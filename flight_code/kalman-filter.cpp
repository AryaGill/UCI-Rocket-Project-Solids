/**
* Implementation of KalmanFilter class.
*
* @author: Dhruv Shah, Hayk Martirosyan
* @date: 07/03/2018
*/

#include <iostream>
#include "kalman-filter.hpp"

KalmanFilter::KalmanFilter(
        const Eigen::MatrixXf& A,
        const Eigen::MatrixXf& B,
        const Eigen::MatrixXf& C,
        const Eigen::MatrixXf& Q,
        const Eigen::MatrixXf& R,
        const Eigen::MatrixXf& P)
    : A(A), B(B), C(C), Q(Q), R(R), P0(P),
      m(C.rows()), n(A.rows()), c(B.cols()), initialized(false),
      I(n, n), x_hat(n)
{
    I.setIdentity();
}

void KalmanFilter::init(const Eigen::VectorXf& x0) {
    x_hat = x0;
    P = P0;
    initialized = true;
}

void KalmanFilter::init() {
    x_hat.setZero();
    P = P0;
    initialized = true;
}

void KalmanFilter::predict(const Eigen::VectorXf& u) {
    if(!initialized) {
        std::cout << "Filter is not initialized! Initializing with trivial state.";
        init();
    }

    x_hat = A * x_hat + B * u;
    P = A * P * A.transpose() + Q;
}

void KalmanFilter::update(const Eigen::VectorXf& y) {
    float S = (C * P * C.transpose())(0,0) + R(0,0);
    K = P * C.transpose() / S;

    x_hat += K * (y - C * x_hat);
    P = (I - K * C) * P;
}

void KalmanFilter::update_dynamics(const Eigen::MatrixXf& A) {
    this->A = A;
}

void KalmanFilter::update_output(const Eigen::MatrixXf& C) {
    this->C = C;
}

void KalmanFilter::update_process_noise(const Eigen::MatrixXf& Q_new) {
    this->Q = Q_new;
}
