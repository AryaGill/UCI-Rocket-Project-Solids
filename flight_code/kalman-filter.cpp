#include "kalman-filter.hpp"
#include <Arduino.h>

// =====================================================
//   Constructor
// =====================================================
KalmanFilter::KalmanFilter(const MatA& A,
                           const MatB& B,
                           const MatC& C,
                           const MatQ& Q,
                           const MatR& R,
                           const MatP& P0)
    : A(A), B(B), C(C), Q(Q), R(R), P0(P0), initialized(false)
{
    I.setIdentity();
    x_hat.setZero();
    P = P0;
}

// Default constructor (initialize to identity)
KalmanFilter::KalmanFilter()
    : initialized(false)
{
    A.setIdentity();
    B.setZero();
    C.setIdentity();
    Q.setIdentity();
    R.setIdentity();
    P0.setIdentity();
    P = P0;
    I.setIdentity();
    x_hat.setZero();
}

// =====================================================
//   Initialization
// =====================================================
void KalmanFilter::init(const VecX& x0) {
    x_hat = x0;
    P = P0;
    initialized = true;
}

void KalmanFilter::init() {
    x_hat.setZero();
    P = P0;
    initialized = true;
}

// =====================================================
//   Predict Step
// =====================================================
void KalmanFilter::predict(const VecU& u) {
    if (!initialized) init();

    // x̂ = A*x + B*u
    x_hat = A * x_hat + B * u;

    // P = A P Aᵀ + Q
    P = A * P * A.transpose() + Q;
}

// =====================================================
//   Update Step
// =====================================================
void KalmanFilter::update(const VecY& y) {
    // Innovation covariance S = C P Cᵀ + R   (3×3)
    Matrix<float,3,3> S = C * P * C.transpose() + R;

    // Kalman gain K = P Cᵀ S⁻¹    (3×3)
    K = P * C.transpose() * S.inverse();

    // Innovation
    VecY innovation = y - C * x_hat;

    // Updated state
    x_hat += K * innovation;

    // Updated covariance
    P = (I - K * C) * P;
}

// =====================================================
//   Matrix update functions
// =====================================================
void KalmanFilter::update_dynamics(const MatA& A_new) {
    A = A_new;
}

void KalmanFilter::update_output(const MatC& C_new) {
    C = C_new;
}

void KalmanFilter::update_process_noise(const MatQ& Q_new) {
    Q = Q_new;
}
