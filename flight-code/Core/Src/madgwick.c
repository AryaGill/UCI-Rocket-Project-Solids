#include "madgwick.h"
#include <math.h>
#include <float.h>

// TODO: Double check these
#define ACCEL_LOW_G   (0.75f * 9.81f)
#define ACCEL_HIGH_G  (2.00f * 9.81f)

static float beta;
static uint32_t prev_time_madgwick = 0;

volatile int g_accel_valid = 0;



static float invSqrt(float x)
{
    if (x <= 0.0f || !isfinite(x)) {
        return 1.0f; // safe default, caller should have guarded
    }
    return 1.0f / sqrtf(x);
}

void Madgwick_Init(Telemetry_t* telemetry, float b)
{
    beta = b;

    float ax = telemetry->lsm_accel_p;
    float ay = telemetry->lsm_accel_y;
    float az = telemetry->lsm_accel_r;

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

    float gx = telemetry->lsm_gyro_p; // rad/s
    float gy = telemetry->lsm_gyro_y;
    float gz = telemetry->lsm_gyro_r;

    float ax = telemetry->lsm_accel_p;
    float ay = telemetry->lsm_accel_y;
    float az = telemetry->lsm_accel_r;

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



//#include "madgwick.h"
//#include <math.h>
//
//#define ACCEL_LOW_G  (0.75f * 9.81f) // 0.75 -> 2 g Fuse accel and mag
//#define ACCEL_HIGH_G (2.00f * 9.81f) // Greater than 2. Gyro only
//
//static float beta;
//static uint32_t prev_time_madgwick = 0;
//
//static float invSqrt(float x){
//	return 1.0f / sqrt(x);
//}
//
//void Madgwick_Init(Telemetry_t* telemetry, float b){
//	beta = b;
//
//	float ax = telemetry->lsm_accel_p;
//	float ay = telemetry->lsm_accel_y;
//	float az = telemetry->lsm_accel_r;
//
//	// Normalize accel
//	float norm = invSqrt(ax*ax + ay*ay + az*az);
//	ax *= norm;
//	ay *= norm;
//	az *= norm;
//
//	 // Compute roll & pitch from gravity
//	float roll  = atan2f(ay, az);
//	float pitch = atan2f(-ax, sqrtf(ay*ay + az*az));
//
//	// Convert to quaternion (yaw = 0)
//	float cr = cosf(roll  * 0.5f);
//	float sr = sinf(roll  * 0.5f);
//	float cp = cosf(pitch * 0.5f);
//	float sp = sinf(pitch * 0.5f);
//
//	telemetry->q0 =  cr * cp;
//	telemetry->q1 =  sr * cp;
//	telemetry->q2 =  cr * sp;
//	telemetry->q3 = -sr * sp;
//
//	prev_time_madgwick = micros();
//}
//
//
//void Madgwick_Update(Telemetry_t* telemetry){
//	// --- dt ---
//	uint32_t cur_time = micros();
//	float dt = (cur_time - prev_time_madgwick) * 1e-6f;
//	prev_time_madgwick = cur_time;
//	if (dt <= 0.0f || dt > 0.05f) return;
//
//	// --- Load state ---
//	float q0 = telemetry->q0;
//	float q1 = telemetry->q1;
//	float q2 = telemetry->q2;
//	float q3 = telemetry->q3;
//
//	float gx = telemetry->lsm_gyro_p;
//	float gy = telemetry->lsm_gyro_y;
//	float gz = telemetry->lsm_gyro_r;
//
//	float ax = telemetry->lsm_accel_p;
//	float ay = telemetry->lsm_accel_y;
//	float az = telemetry->lsm_accel_r;
//
//	float mx = telemetry->mag_p;
//	float my = telemetry->mag_y;
//	float mz = telemetry->mag_r;
//
//	// --- Gyro quaternion derivative (always applied) ---
//	float qDot0 = 0.5f * (-q1*gx - q2*gy - q3*gz);
//	float qDot1 = 0.5f * ( q0*gx + q2*gz - q3*gy);
//	float qDot2 = 0.5f * ( q0*gy - q1*gz + q3*gx);
//	float qDot3 = 0.5f * ( q0*gz + q1*gy - q2*gx);
//
//
//	// ACCEL GATE
//	float accel_mag = sqrtf(ax*ax + ay*ay + az*az);
//    int accel_valid = (accel_mag > ACCEL_LOW_G && accel_mag < ACCEL_HIGH_G);
//
//	// --- Mag validity ---
//	float mag_norm = mx*mx + my*my + mz*mz;
//	int mag_valid = (mag_norm > 1e-6f) && accel_valid;
//
//	if (accel_valid)
//	{
//		float recipNorm;
//		float s0, s1, s2, s3;
//
//		// Normalize accel
//		recipNorm = invSqrt(ax*ax + ay*ay + az*az);
//		ax *= recipNorm;
//		ay *= recipNorm;
//		az *= recipNorm;
//
//		if (mag_valid)
//		{
//			// Normalize mag
//			recipNorm = invSqrt(mag_norm);
//			mx *= recipNorm;
//			my *= recipNorm;
//			mz *= recipNorm;
//
//			// Reference direction of Earth's magnetic field
//			float _2q0mx = 2.0f*q0*mx;
//			float _2q0my = 2.0f*q0*my;
//			float _2q0mz = 2.0f*q0*mz;
//			float _2q1mx = 2.0f*q1*mx;
//			float _2q0 = 2.0f*q0, _2q1 = 2.0f*q1;
//			float _2q2 = 2.0f*q2, _2q3 = 2.0f*q3;
//			float q0q0=q0*q0, q1q1=q1*q1, q2q2=q2*q2, q3q3=q3*q3;
//
//			float hx = mx*q0q0 - _2q0my*q3 + _2q0mz*q2
//					 + mx*q1q1 + _2q1*my*q2 + _2q1*mz*q3
//					 - mx*q2q2 - mx*q3q3;
//
//			float hy = _2q0*mx*q3 + my*q0q0 - _2q0mz*q1
//					 + _2q1*mx*q2 - my*q1q1 + my*q2q2
//					 + _2q2*mz*q3 - my*q3q3;
//
//			float _2bx = sqrtf(hx*hx + hy*hy);
//			float _2bz = -_2q0*mx*q2 + _2q0*my*q1 + mz*q0q0
//					   + _2q1*mx*q3 - mz*q1q1
//					   + _2q2*my*q3 - mz*q2q2 + mz*q3q3;
//
//			// Gradient descent — accel + mag
//			s0 = -_2q2*(2*(q1*q3 - q0*q2) - ax)
//			   + _2q1*(2*(q0*q1 + q2*q3) - ay)
//			   - _2bz*q2*(_2bx*(0.5f-q2q2-q3q3) + _2bz*(q1*q3-q0*q2) - mx)
//			   + (-_2bx*q3+_2bz*q1)*(_2bx*(q1*q2-q0*q3) + _2bz*(q0*q1+q2*q3) - my)
//			   + _2bx*q2*(_2bx*(q0*q2+q1*q3) + _2bz*(0.5f-q1q1-q2q2) - mz);
//
//			s1 = _2q3*(2*(q1*q3 - q0*q2) - ax)
//			   + _2q0*(2*(q0*q1 + q2*q3) - ay)
//			   - 4*q1*(1 - 2*(q1q1+q2q2) - az)
//			   + _2bz*q3*(_2bx*(0.5f-q2q2-q3q3) + _2bz*(q1*q3-q0*q2) - mx)
//			   + (_2bx*q2+_2bz*q0)*(_2bx*(q1*q2-q0*q3) + _2bz*(q0*q1+q2*q3) - my)
//			   + (_2bx*q3-4*_2bz*q1)*(_2bx*(q0*q2+q1*q3) + _2bz*(0.5f-q1q1-q2q2) - mz);
//
//			s2 = -_2q0*(2*(q1*q3 - q0*q2) - ax)
//			   + _2q3*(2*(q0*q1 + q2*q3) - ay)
//			   - 4*q2*(1 - 2*(q1q1+q2q2) - az)
//			   + (-4*_2bx*q2-_2bz*q0)*(_2bx*(0.5f-q2q2-q3q3) + _2bz*(q1*q3-q0*q2) - mx)
//			   + (_2bx*q1+_2bz*q3)*(_2bx*(q1*q2-q0*q3) + _2bz*(q0*q1+q2*q3) - my)
//			   + (_2bx*q0-4*_2bz*q2)*(_2bx*(q0*q2+q1*q3) + _2bz*(0.5f-q1q1-q2q2) - mz);
//
//			s3 = _2q1*(2*(q1*q3 - q0*q2) - ax)
//			   + _2q2*(2*(q0*q1 + q2*q3) - ay)
//			   + (-4*_2bx*q3+_2bz*q1)*(_2bx*(0.5f-q2q2-q3q3) + _2bz*(q1*q3-q0*q2) - mx)
//			   + (-_2bx*q0+_2bz*q2)*(_2bx*(q1*q2-q0*q3) + _2bz*(q0*q1+q2*q3) - my)
//			   + _2bx*q1*(_2bx*(q0*q2+q1*q3) + _2bz*(0.5f-q1q1-q2q2) - mz);
//		}
//		else
//		{
//			// Gradient descent — accel only (no mag)
//			float q0q0=q0*q0, q1q1=q1*q1, q2q2=q2*q2, q3q3=q3*q3;
//
//			s0 =  4*q0*q2*q2 + 2*q2*ax + 4*q0*q1*q1 - 2*q1*ay;
//			s1 =  4*q1*q3*q3 - 2*q3*ax + 4*q0*q0*q1
//				- 2*q0*ay - 4*q1 + 8*q1*q1q1
//				+ 8*q1*q2q2 + 4*q1*az;
//			s2 =  4*q0*q0*q2 + 2*q0*ax + 4*q2*q3*q3
//				- 2*q3*ay - 4*q2 + 8*q2*q1q1
//				+ 8*q2*q2q2 + 4*q2*az;
//			s3 =  4*q1*q1*q3 - 2*q1*ax + 4*q2*q2*q3 - 2*q2*ay;
//		}
//
//		// Normalize gradient
//		float norm = s0*s0 + s1*s1 + s2*s2 + s3*s3;
//		if (norm > 0.0f)
//		{
//			recipNorm = invSqrt(norm);
//			s0 *= recipNorm;
//			s1 *= recipNorm;
//			s2 *= recipNorm;
//			s3 *= recipNorm;
//		}
//
//		// Apply correction
//		qDot0 -= beta * s0;
//		qDot1 -= beta * s1;
//		qDot2 -= beta * s2;
//		qDot3 -= beta * s3;
//	}
//	// else: high-G burn → gyro only, no correction applied
//
//	// --- Integrate ---
//	q0 += qDot0 * dt;
//	q1 += qDot1 * dt;
//	q2 += qDot2 * dt;
//	q3 += qDot3 * dt;
//
//	// --- Renormalize ---
//	float recipNorm = invSqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
//	telemetry->q0 = q0 * recipNorm;
//	telemetry->q1 = q1 * recipNorm;
//	telemetry->q2 = q2 * recipNorm;
//	telemetry->q3 = q3 * recipNorm;
//}
//
//void Madgwick_GetEuler(Telemetry_t* telemetry)
//{
//    float q0 = telemetry->q0;
//    float q1 = telemetry->q1;
//    float q2 = telemetry->q2;
//    float q3 = telemetry->q3;
//
//    // Roll (rotation about body X)
//    float sinr_cosp = 2.0f * (q0*q1 + q2*q3);
//    float cosr_cosp = 1.0f - 2.0f * (q1*q1 + q2*q2);
//    telemetry->roll = atan2f(sinr_cosp, cosr_cosp);
//
//    // Pitch (rotation about body Y)
//    float sinp = 2.0f * (q0*q2 - q3*q1);
//    if (fabsf(sinp) >= 1.0f)
//        telemetry->pitch = copysignf(M_PI / 2.0f, sinp); // clamp at ±90°
//    else
//        telemetry->pitch = asinf(sinp);
//
//    // Yaw (rotation about body Z)
//    float siny_cosp = 2.0f * (q0*q3 + q1*q2);
//    float cosy_cosp = 1.0f - 2.0f * (q2*q2 + q3*q3);
//    telemetry->yaw = atan2f(siny_cosp, cosy_cosp);
//
//    // Convert to degrees
//    telemetry->roll  *= (180.0f / M_PI);
//    telemetry->pitch *= (180.0f / M_PI);
//    telemetry->yaw   *= (180.0f / M_PI);
//}
//
//
////#include "madgwick.h"
////#include <math.h>
////
////float beta;
////uint32_t prev_time_madgwick = 0;
////
////static float invSqrt(float x)
////{
////    return 1.0f / sqrtf(x);
////}
////
////
////void Madgwick_Init(Telemetry_t* telemetry, float b)
////{
////    beta = b;
////
////    float ax = telemetry->lsm_accel_p;
////	float ay = telemetry->lsm_accel_y;
////	float az = telemetry->lsm_accel_r;
////
////    // Normalize accel
////    float norm = sqrtf(ax*ax + ay*ay + az*az);
////    ax /= norm;
////    ay /= norm;
////    az /= norm;
////
////    // Compute roll & pitch
////    float roll  = atan2f(ay, az);
////    float pitch = atan2f(-ax, sqrtf(ay*ay + az*az));
////
////    // Convert to quaternion (yaw = 0)
////    float cr = cosf(roll * 0.5f);
////    float sr = sinf(roll * 0.5f);
////    float cp = cosf(pitch * 0.5f);
////    float sp = sinf(pitch * 0.5f);
////
////    telemetry->q0 = cr * cp;
////    telemetry->q1 = sr * cp;
////    telemetry->q2 = cr * sp;
////    telemetry->q3 = sr * sp;
////
////    prev_time_madgwick = micros();
////}
////
////
/////* ================= IMU UPDATE ================= */
////float dt;
////void Madgwick_UpdateIMU(Telemetry_t* telemetry)
////{
////	// get dt (time in ms since last update)
////	uint64_t cur_time = micros();
////	dt = (cur_time - prev_time_madgwick) * 1e-6f;
////	prev_time_madgwick = cur_time;
////	if (dt <= 0.0f || dt > 0.05f) return;
////
////	float q0 = telemetry->q0;
////	float q1 = telemetry->q1;
////	float q2 = telemetry->q2;
////	float q3 = telemetry->q3;
////
////    float ax = telemetry->lsm_accel_p;
////	float ay = telemetry->lsm_accel_y;
////	float az = telemetry->lsm_accel_r;
////
////	float gx = telemetry->lsm_gyro_p;
////	float gy = telemetry->lsm_gyro_y;
////	float gz = telemetry->lsm_gyro_r;
////
////    float recipNorm;
////    float s0, s1, s2, s3;
////    float qDot1, qDot2, qDot3, qDot4;
////
////
////    /* Gyro quaternion derivative */
////    qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
////    qDot2 = 0.5f * ( q0 * gx + q2 * gz - q3 * gy);
////    qDot3 = 0.5f * ( q0 * gy - q1 * gz + q3 * gx);
////    qDot4 = 0.5f * ( q0 * gz + q1 * gy - q2 * gx);
////
////
////    if (!(ax == 0.0f && ay == 0.0f && az == 0.0f))
////    {
////        recipNorm = invSqrt(ax*ax + ay*ay + az*az);
////        ax *= recipNorm;
////        ay *= recipNorm;
////        az *= recipNorm;
////
////
////        s0 = 4.0f*q0*q2*q2 + 2.0f*q2*ax + 4.0f*q0*q1*q1 - 2.0f*q1*ay;
////        s1 = 4.0f*q1*q3*q3 - 2.0f*q3*ax + 4.0f*q0*q0*q1
////           - 2.0f*q0*ay - 4.0f*q1 + 8.0f*q1*q1*q1
////           + 8.0f*q1*q2*q2 + 4.0f*q1*az;
////
////        s2 = 4.0f*q0*q0*q2 + 2.0f*q0*ax + 4.0f*q2*q3*q3
////           - 2.0f*q3*ay - 4.0f*q2 + 8.0f*q2*q1*q1
////           + 8.0f*q2*q2*q2 + 4.0f*q2*az;
////
////        s3 = 4.0f*q1*q1*q3 - 2.0f*q1*ax + 4.0f*q2*q2*q3 - 2.0f*q2*ay;
////
////        float norm = s0*s0 + s1*s1 + s2*s2 + s3*s3;
////		if (norm > 0.0f) {
////			recipNorm = invSqrt(norm);
////			s0 *= recipNorm;
////			s1 *= recipNorm;
////			s2 *= recipNorm;
////			s3 *= recipNorm;
////		}
////
////        qDot1 -= beta * s0;
////        qDot2 -= beta * s1;
////        qDot3 -= beta * s2;
////        qDot4 -= beta * s3;
////    }
////
////
////    /* Integrate */
////    q0 += qDot1 * dt;
////    q1 += qDot2 * dt;
////    q2 += qDot3 * dt;
////    q3 += qDot4 * dt;
////
////    recipNorm = invSqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
////
////    telemetry->q0 = q0 * recipNorm;
////    telemetry->q1 = q1 * recipNorm;
////    telemetry->q2 = q2 * recipNorm;
////    telemetry->q3 = q3 * recipNorm;
////}
////
////
////
/////* ================= FULL AHRS ================= */
////
////void Madgwick_Update(Telemetry_t* telemetry)
////{
//////	Madgwick_UpdateIMU(telemetry);
//////	return;
////
////	float mx = telemetry->mag_p;
////	float my = telemetry->mag_y;
////	float mz = telemetry->mag_r;
////
////	float mag_norm = mx*mx + my*my + mz*mz;
////
////	if (mag_norm < 1e-6f) {
////	    Madgwick_UpdateIMU(telemetry);
////	    return;
////	}
////
////    // get dt (time in ms since last update)
////    uint64_t cur_time = micros();
////	float dt = (cur_time - prev_time_madgwick) * 1e-6f;
////	prev_time_madgwick = cur_time;
////	if (dt <= 0.0f || dt > 0.05f) return;
////
////    float q0 = telemetry->q0;
////    float q1 = telemetry->q1;
////    float q2 = telemetry->q2;
////    float q3 = telemetry->q3;
////
////    float ax = telemetry->lsm_accel_p;
////    float ay = telemetry->lsm_accel_y;
////    float az = telemetry->lsm_accel_r;
////
////    float gx = telemetry->lsm_gyro_p;
////	float gy = telemetry->lsm_gyro_y;
////	float gz = telemetry->lsm_gyro_r;
////
////	float recipNorm;
////	float s0, s1, s2, s3;
////	float qDot1, qDot2, qDot3, qDot4;
////
////	float hx, hy;
////	float _2bx, _2bz;
////
////	float _2q0mx, _2q0my, _2q0mz, _2q1mx;
////	float _2q0 = 2.0f * q0;
////	float _2q1 = 2.0f * q1;
////	float _2q2 = 2.0f * q2;
////	float _2q3 = 2.0f * q3;
////
////	float q0q0 = q0*q0;
////	float q1q1 = q1*q1;
////	float q2q2 = q2*q2;
////	float q3q3 = q3*q3;
////
////	/* normalize accel */
////	float norm = ax*ax + ay*ay + az*az;
////	if (norm == 0.0f) return;
////
////	recipNorm = invSqrt(norm);
////	ax *= recipNorm;
////	ay *= recipNorm;
////	az *= recipNorm;
////
////	/* normalize mag */
////	recipNorm = invSqrt(mx*mx + my*my + mz*mz);
////	mx *= recipNorm;
////	my *= recipNorm;
////	mz *= recipNorm;
////
////	/* reference direction of Earth's magnetic field */
////	_2q0mx = 2.0f*q0*mx;
////	_2q0my = 2.0f*q0*my;
////	_2q0mz = 2.0f*q0*mz;
////	_2q1mx = 2.0f*q1*mx;
////
////	hx = mx*q0q0 - _2q0my*q3 + _2q0mz*q2 +
////		 mx*q1q1 + _2q1*my*q2 + _2q1*mz*q3 -
////		 mx*q2q2 - mx*q3q3;
////
////	hy = _2q0*mx*q3 + my*q0q0 - _2q0mz*q1 +
////		 _2q1*mx*q2 - my*q1q1 + my*q2q2 +
////		 _2q2*mz*q3 - my*q3q3;
////
////	_2bx = sqrtf(hx*hx + hy*hy);
////	_2bz = -_2q0*mx*q2 + _2q0*my*q1 + mz*q0q0 +
////		   _2q1*mx*q3 - mz*q1q1 +
////		   _2q2*my*q3 - mz*q2q2 + mz*q3q3;
////
////
////    /* Gyro derivative */
////    qDot1 = 0.5f * (-q1*gx - q2*gy - q3*gz);
////    qDot2 = 0.5f * ( q0*gx + q2*gz - q3*gy);
////    qDot3 = 0.5f * ( q0*gy - q1*gz + q3*gx);
////    qDot4 = 0.5f * ( q0*gz + q1*gy - q2*gx);
////
////
////    /* Gradient descent step (trimmed but correct) */
//////    s0 = -2.0f*q2*(2*(q1*q3 - q0*q2) - ax)
//////       + 2.0f*q1*(2*(q0*q1 + q2*q3) - ay)
//////       - _2bz*q2*(_2bx*(0.5f - q2q2 - q3q3) + _2bz*(q1*q3 - q0*q2) - mx);
//////
//////    s1 =  2.0f*q3*(2*(q1*q3 - q0*q2) - ax)
//////       + 2.0f*q0*(2*(q0*q1 + q2*q3) - ay)
//////       - 4.0f*q1*(1 - 2*(q1q1 + q2q2) - az);
//////
//////    s2 = -2.0f*q0*(2*(q1*q3 - q0*q2) - ax)
//////       + 2.0f*q3*(2*(q0*q1 + q2*q3) - ay)
//////       - 4.0f*q2*(1 - 2*(q1q1 + q2q2) - az);
//////
//////    s3 =  2.0f*q1*(2*(q1*q3 - q0*q2) - ax)
//////       + 2.0f*q2*(2*(q0*q1 + q2*q3) - ay);
////
////    s0 = -_2q2*(2*(q1*q3 - q0*q2) - ax)
////		 + _2q1*(2*(q0*q1 + q2*q3) - ay)
////		 - _2bz*q2*(_2bx*(0.5f - q2q2 - q3q3) + _2bz*(q1*q3 - q0*q2) - mx)
////		 + (-_2bx*q3 + _2bz*q1)*(_2bx*(q1*q2 - q0*q3) + _2bz*(q0*q1 + q2*q3) - my)
////		 + _2bx*q2*(_2bx*(q0*q2 + q1*q3) + _2bz*(0.5f - q1q1 - q2q2) - mz);
////
////	s1 = _2q3*(2*(q1*q3 - q0*q2) - ax)
////		 + _2q0*(2*(q0*q1 + q2*q3) - ay)
////		 - 4*q1*(1 - 2*(q1q1 + q2q2) - az)
////		 + _2bz*q3*(_2bx*(0.5f - q2q2 - q3q3) + _2bz*(q1*q3 - q0*q2) - mx)
////		 + (_2bx*q2 + _2bz*q0)*(_2bx*(q1*q2 - q0*q3) + _2bz*(q0*q1 + q2*q3) - my)
////		 + (_2bx*q3 - 4*_2bz*q1)*(_2bx*(q0*q2 + q1*q3) + _2bz*(0.5f - q1q1 - q2q2) - mz);
////
////	s2 = -_2q0*(2*(q1*q3 - q0*q2) - ax)
////		 + _2q3*(2*(q0*q1 + q2*q3) - ay)
////		 - 4*q2*(1 - 2*(q1q1 + q2q2) - az)
////		 + (-4*_2bx*q2 - _2bz*q0)*(_2bx*(0.5f - q2q2 - q3q3) + _2bz*(q1*q3 - q0*q2) - mx)
////		 + (_2bx*q1 + _2bz*q3)*(_2bx*(q1*q2 - q0*q3) + _2bz*(q0*q1 + q2*q3) - my)
////		 + (_2bx*q0 - 4*_2bz*q2)*(_2bx*(q0*q2 + q1*q3) + _2bz*(0.5f - q1q1 - q2q2) - mz);
////
////	s3 = _2q1*(2*(q1*q3 - q0*q2) - ax)
////		 + _2q2*(2*(q0*q1 + q2*q3) - ay)
////		 + (-4*_2bx*q3 + _2bz*q1)*(_2bx*(0.5f - q2q2 - q3q3) + _2bz*(q1*q3 - q0*q2) - mx)
////		 + (-_2bx*q0 + _2bz*q2)*(_2bx*(q1*q2 - q0*q3) + _2bz*(q0*q1 + q2*q3) - my)
////		 + _2bx*q1*(_2bx*(q0*q2 + q1*q3) + _2bz*(0.5f - q1q1 - q2q2) - mz);
////
////
////
////
////    norm = s0*s0 + s1*s1 + s2*s2 + s3*s3;
////    if (norm > 0.0f) {
////        recipNorm = invSqrt(norm);
////        s0 *= recipNorm;
////        s1 *= recipNorm;
////        s2 *= recipNorm;
////        s3 *= recipNorm;
////    }
////
////    qDot1 -= beta * s0;
////    qDot2 -= beta * s1;
////    qDot3 -= beta * s2;
////    qDot4 -= beta * s3;
////
////
////    /* Integrate */
////    q0 += qDot1 * dt;
////    q1 += qDot2 * dt;
////    q2 += qDot3 * dt;
////    q3 += qDot4 * dt;
////
////    recipNorm = invSqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
////
////    telemetry->q0 = q0 * recipNorm;
////    telemetry->q1 = q1 * recipNorm;
////    telemetry->q2 = q2 * recipNorm;
////    telemetry->q3 = q3 * recipNorm;
////}
