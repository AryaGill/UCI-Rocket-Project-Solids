#include "bias.h"

void Bias_Init(Bias_t *bias)
{
    bias->adxl_accel_r_bias = 0.0f;
    bias->adxl_accel_p_bias = 0.0f;
    bias->adxl_accel_y_bias = 0.0f;

    bias->lsm_accel_r_bias = 0.0f;
    bias->lsm_accel_p_bias = 0.0f;
    bias->lsm_accel_y_bias = 0.0f;

    bias->bias_count = 0;
}

void Bias_Calculate(Bias_t *bias, Telemetry_t *t)
{
    bias->bias_count += 1.0f;

    float n = bias->bias_count;

    bias->adxl_accel_r_bias += ((t->adxl_accel_r - 9.81) - bias->adxl_accel_r_bias) / n;
    bias->adxl_accel_p_bias += (t->adxl_accel_p - bias->adxl_accel_p_bias) / n;
    bias->adxl_accel_y_bias += (t->adxl_accel_y - bias->adxl_accel_y_bias) / n;

    bias->lsm_accel_r_bias += ((t->lsm_accel_r - 9.81) - bias->lsm_accel_r_bias) / n;
    bias->lsm_accel_p_bias += (t->lsm_accel_p - bias->lsm_accel_p_bias) / n;
    bias->lsm_accel_y_bias += (t->lsm_accel_y - bias->lsm_accel_y_bias) / n;
}

void Apply_Bias(Bias_t *bias, Telemetry_t *t)
{
    t->adxl_accel_r -= bias->adxl_accel_r_bias;
    t->adxl_accel_p -= bias->adxl_accel_p_bias;
    t->adxl_accel_y -= bias->adxl_accel_y_bias;

    t->lsm_accel_r -= bias->lsm_accel_r_bias;
    t->lsm_accel_p -= bias->lsm_accel_p_bias;
    t->lsm_accel_y -= bias->lsm_accel_y_bias;
}
