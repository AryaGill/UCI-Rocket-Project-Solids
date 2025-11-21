#include <ArduinoEigen.h>
#include <ArduinoEigenDense.h>
#include <ArduinoEigenSparse.h>

#pragma once

// ======================================================================
// Eigen configuration for embedded / microcontroller builds
// ======================================================================

// Disable Eigen features that cause issues on ARM/Teensy
#include <ArduinoEigenDense.h>

using namespace Eigen;

// ======================================================================
// KalmanFilter class definition
// ======================================================================

class KalmanFilter {
public:
    KalmanFilter(
        const Eigen::MatrixXf& A,
        const Eigen::MatrixXf& B,
        const Eigen::MatrixXf& C,
        const Eigen::MatrixXf& Q,
        const Eigen::MatrixXf& R,
        const Eigen::MatrixXf& P);

    KalmanFilter();

    void init(const Eigen::VectorXf& x0);
    void init();

    void predict(const Eigen::VectorXf& u);
    void update(const Eigen::VectorXf& y);

    void update_dynamics(const Eigen::MatrixXf& A);
    void update_output(const Eigen::MatrixXf& C);
    void update_process_noise(const Eigen::MatrixXf& Q_new);

    Eigen::VectorXf state() const { return x_hat; }

private:
    Eigen::MatrixXf A, B, C, Q, R, P, K, P0;
    int m, n, c;
    bool initialized;

    Eigen::MatrixXf I;
    Eigen::VectorXf x_hat;
};
