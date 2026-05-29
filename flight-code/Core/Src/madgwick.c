#include "madgwick.h"
#include <math.h>
#include <float.h>

// TODO: Double check these
#define ACCEL_LOW_G   (0.75f * 9.81f)
#define ACCEL_HIGH_G  (2.00f * 9.81f)

static float beta;
static uint32_t prev_time_madgwick = 0;

volatile int g_accel_valid = 0;


//compute inverse square root of a value
static float invSqrt(float x)
{
    if (x <= 0.0f || !isfinite(x)) {
        return 1.0f; // safe default, caller should have guarded
    }
    return 1.0f / sqrtf(x);
}

/**
 * Initializes the Madgwick filter.
 *
 * telemetry: flight telemetry struct
 * b: filter gain value
 *
 * Sets the initial quaternion from accelerometer data.
 */
void Madgwick_Init(Telemetry_t* telemetry, float b)
{
    beta = b;

    float ax = telemetry->bmx_accel_p;
    float ay = telemetry->bmx_accel_y;
    float az = telemetry->bmx_accel_r;

    // guard accel norm
    float norm = ax*ax + ay*ay + az*az;
    if (norm <= 0.0f || !isfinite(norm)) {
        // fall back to identity quaternion
        telemetry->q0 = 1.0f;
        telemetry->q1 = 0.0f;
        telemetry->q2 = 0.0f;
        telemetry->q3 = 0.0f;
        prev_time_madgwick = micros();
        return;
    }

    float recipNorm = invSqrt(norm);
    ax *= recipNorm;
    ay *= recipNorm;
    az *= recipNorm;

    // Compute roll & pitch from gravity
    float roll  = atan2f(ay, az);
    float pitch = atan2f(-ax, sqrtf(ay*ay + az*az));

    // yaw = 0
    float cr = cosf(roll  * 0.5f);
    float sr = sinf(roll  * 0.5f);
    float cp = cosf(pitch * 0.5f);
    float sp = sinf(pitch * 0.5f);

    telemetry->q0 =  cr * cp;
    telemetry->q1 =  sr * cp;
    telemetry->q2 =  cr * sp;
    telemetry->q3 = -sr * sp;

    prev_time_madgwick = micros();
}

/**
 * Updates the Madgwick filter.
 *
 * telemetry: flight telemetry struct
 *
 * Uses gyroscope and accelerometer data to update the orientation quaternion.
 */
void Madgwick_Update(Telemetry_t* telemetry)
{
    // --- dt ---
    uint32_t cur_time = micros();
    float dt = (cur_time - prev_time_madgwick) * 1e-6f;
    prev_time_madgwick = cur_time;

    if (!isfinite(dt) || dt <= 0.0f || dt > 0.05f) {
        return;
    }

    // --- Load state ---
    float q0 = telemetry->q0;
    float q1 = telemetry->q1;
    float q2 = telemetry->q2;
    float q3 = telemetry->q3;

    float gx = telemetry->bmx_gyro_p; // rad/s
    float gy = telemetry->bmx_gyro_y;
    float gz = telemetry->bmx_gyro_r;

    float ax = telemetry->bmx_accel_p;
    float ay = telemetry->bmx_accel_y;
    float az = telemetry->bmx_accel_r;

    // basic sanity on inputs
    if (!isfinite(q0) || !isfinite(q1) || !isfinite(q2) || !isfinite(q3) ||
        !isfinite(gx) || !isfinite(gy) || !isfinite(gz) ||
        !isfinite(ax) || !isfinite(ay) || !isfinite(az)) {
        // reset to identity if things blew up
        telemetry->q0 = 1.0f;
        telemetry->q1 = 0.0f;
        telemetry->q2 = 0.0f;
        telemetry->q3 = 0.0f;
        return;
    }

    // --- Gyro quaternion derivative (always applied) ---
    float qDot0 = 0.5f * (-q1*gx - q2*gy - q3*gz);
    float qDot1 = 0.5f * ( q0*gx + q2*gz - q3*gy);
    float qDot2 = 0.5f * ( q0*gy - q1*gz + q3*gx);
    float qDot3 = 0.5f * ( q0*gz + q1*gy - q2*gx);

    // --- Optional accel correction (IMU-only Madgwick) ---
    float accel_mag = sqrtf(ax*ax + ay*ay + az*az);
    int accel_valid = (accel_mag > ACCEL_LOW_G && accel_mag < ACCEL_HIGH_G);
    g_accel_valid = accel_valid;

    if (accel_valid) {
        // normalize accel
        float recipNorm = invSqrt(ax*ax + ay*ay + az*az);
        ax *= recipNorm;
        ay *= recipNorm;
        az *= recipNorm;

        // reference direction of gravity from quaternion
        float _2q0 = 2.0f * q0;
        float _2q1 = 2.0f * q1;
        float _2q2 = 2.0f * q2;
        float _2q3 = 2.0f * q3;
        float _4q0 = 4.0f * q0;
        float _4q1 = 4.0f * q1;
        float _4q2 = 4.0f * q2;
        float _8q1 = 8.0f * q1;
        float _8q2 = 8.0f * q2;
        float q0q0 = q0 * q0;
        float q1q1 = q1 * q1;
        float q2q2 = q2 * q2;
        float q3q3 = q3 * q3;

        // gradient of objective function (standard Madgwick IMU form)
        float s0 = _4q0*q2q2 + _2q2*ax + _4q0*q1q1 - _2q1*ay;
        float s1 = _4q1*q3q3 - _2q3*ax + 4.0f*q0q0*q1 - _2q0*ay
                 - _4q1 + _8q1*q1q1 + _8q1*q2q2 + _4q1*az;
        float s2 = 4.0f*q0q0*q2 + _2q0*ax + _4q2*q3q3 - _2q3*ay
                 - _4q2 + _8q2*q1q1 + _8q2*q2q2 + _4q2*az;
        float s3 = 4.0f*q1q1*q3 - _2q1*ax + 4.0f*q2q2*q3 - _2q2*ay;

        float snorm = s0*s0 + s1*s1 + s2*s2 + s3*s3;
        if (snorm > 0.0f && isfinite(snorm)) {
            recipNorm = invSqrt(snorm);
            s0 *= recipNorm;
            s1 *= recipNorm;
            s2 *= recipNorm;
            s3 *= recipNorm;

            // Apply correction
            qDot0 -= beta * s0;
            qDot1 -= beta * s1;
            qDot2 -= beta * s2;
            qDot3 -= beta * s3;
        }
    }
    // else: high-G → gyro-only

    // --- Integrate ---
    q0 += qDot0 * dt;
    q1 += qDot1 * dt;
    q2 += qDot2 * dt;
    q3 += qDot3 * dt;

    // NaN guard before renormalization
    if (!isfinite(q0) || !isfinite(q1) || !isfinite(q2) || !isfinite(q3)) {
        telemetry->q0 = 1.0f;
        telemetry->q1 = 0.0f;
        telemetry->q2 = 0.0f;
        telemetry->q3 = 0.0f;
        return;
    }

    // --- Renormalize ---
    float qnorm2 = q0*q0 + q1*q1 + q2*q2 + q3*q3;
    if (qnorm2 <= 0.0f || !isfinite(qnorm2)) {
        telemetry->q0 = 1.0f;
        telemetry->q1 = 0.0f;
        telemetry->q2 = 0.0f;
        telemetry->q3 = 0.0f;
        return;
    }

    float recipNorm = invSqrt(qnorm2);
    telemetry->q0 = q0 * recipNorm;
    telemetry->q1 = q1 * recipNorm;
    telemetry->q2 = q2 * recipNorm;
    telemetry->q3 = q3 * recipNorm;
}

//converts current quaternion into euler angles-> roll pitch yaw
void Madgwick_GetEuler(Telemetry_t* telemetry)
{
    float q0 = telemetry->q0;
    float q1 = telemetry->q1;
    float q2 = telemetry->q2;
    float q3 = telemetry->q3;

    if (!isfinite(q0) || !isfinite(q1) || !isfinite(q2) || !isfinite(q3)) {
        telemetry->roll  = 0.0f;
        telemetry->pitch = 0.0f;
        telemetry->yaw   = 0.0f;
        return;
    }

    // Roll (rotation about body X)
    float sinr_cosp = 2.0f * (q0*q1 + q2*q3);
    float cosr_cosp = 1.0f - 2.0f * (q1*q1 + q2*q2);
    telemetry->pitch = atan2f(sinr_cosp, cosr_cosp);

    // Pitch (rotation about body Y)
    float sinp = 2.0f * (q0*q2 - q3*q1);
    if (fabsf(sinp) >= 1.0f)
        telemetry->yaw = copysignf(M_PI / 2.0f, sinp);
    else
        telemetry->yaw = asinf(sinp);

    // Yaw (rotation about body Z)
    float siny_cosp = 2.0f * (q0*q3 + q1*q2);
    float cosy_cosp = 1.0f - 2.0f * (q2*q2 + q3*q3);
    telemetry->roll = atan2f(siny_cosp, cosy_cosp);

    // Convert to degrees
    const float RAD2DEG = 180.0f / (float)M_PI;
    telemetry->roll  *= RAD2DEG;
    telemetry->pitch *= RAD2DEG;
    telemetry->yaw   *= RAD2DEG;
}
