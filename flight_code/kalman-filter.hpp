#pragma once
#include <ArduinoEigen.h>
#include <ArduinoEigenDense.h>

using namespace Eigen;

class KalmanFilter {
public:
    // -------- Fixed sizes --------
    static constexpr int n = 3;   // state dimension
    static constexpr int m = 3;   // measurement dimension
    static constexpr int c = 1;   // control dimension

    using MatA = Matrix<float, n, n>;       // 3×3
    using MatB = Matrix<float, n, c>;       // 3×1
    using MatC = Matrix<float, m, n>;       // 3×3
    using MatP = Matrix<float, n, n>;       // 3×3
    using MatQ = Matrix<float, n, n>;       // 3×3
    using MatR = Matrix<float, m, m>;       // 3×3
    using MatK = Matrix<float, n, m>;       // 3×3
    using VecX = Matrix<float, n, 1>;       // 3×1
    using VecU = Matrix<float, c, 1>;       // 1×1
    using VecY = Matrix<float, m, 1>;       // 3×1

    // -------- Constructor --------
    KalmanFilter(const MatA& A,
                 const MatB& B,
                 const MatC& C,
                 const MatQ& Q,
                 const MatR& R,
                 const MatP& P0);

    // Default constructor
    KalmanFilter();

    // -------- Initialization --------
    void init(const VecX& x0);
    void init();

    // -------- Predict + Update --------
    void predict(const VecU& u);
    void update(const VecY& y);

    void update_dynamics(const MatA& A_new);
    void update_output(const MatC& C_new);
    void update_process_noise(const MatQ& Q_new);

    VecX state() const { return x_hat; }

private:
    MatA A;
    MatB B;
    MatC C;
    MatQ Q;
    MatR R;
    MatP P;
    MatP P0;
    MatP I;
    MatK K;

    VecX x_hat;

    bool initialized;
};
