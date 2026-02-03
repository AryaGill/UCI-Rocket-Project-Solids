#include "madgwick.h"
#include <math.h>

float beta;

float *Q0;
float *Q1;
float *Q2;
float *Q3;

static float invSqrt(float x)
{
    return 1.0f / sqrtf(x);
}


void Madgwick_Init(float *q0, float *q1, float *q2, float *q3, float b)
{
    beta = b;

    Q0 = q0;
    Q1 = q1;
    Q2 = q2;
    Q3 = q3;

    *Q0 = 1.0f;
    *Q1 = 0.0f;
    *Q2 = 0.0f;
    *Q3 = 0.0f;

}


/* ================= IMU UPDATE ================= */

void Madgwick_UpdateIMU(float gx, float gy, float gz,
                        float ax, float ay, float az,
                        float dt)
{
    float q0 = *Q0;
    float q1 = *Q1;
    float q2 = *Q2;
    float q3 = *Q3;

    float recipNorm;
    float s0, s1, s2, s3;
    float qDot1, qDot2, qDot3, qDot4;


    /* Gyro quaternion derivative */
    qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
    qDot2 = 0.5f * ( q0 * gx + q2 * gz - q3 * gy);
    qDot3 = 0.5f * ( q0 * gy - q1 * gz + q3 * gx);
    qDot4 = 0.5f * ( q0 * gz + q1 * gy - q2 * gx);


    if (!(ax == 0.0f && ay == 0.0f && az == 0.0f))
    {
        recipNorm = invSqrt(ax*ax + ay*ay + az*az);
        ax *= recipNorm;
        ay *= recipNorm;
        az *= recipNorm;


        s0 = 4.0f*q0*q2*q2 + 2.0f*q2*ax + 4.0f*q0*q1*q1 - 2.0f*q1*ay;
        s1 = 4.0f*q1*q3*q3 - 2.0f*q3*ax + 4.0f*q0*q0*q1
           - 2.0f*q0*ay - 4.0f*q1 + 8.0f*q1*q1*q1
           + 8.0f*q1*q2*q2 + 4.0f*q1*az;

        s2 = 4.0f*q0*q0*q2 + 2.0f*q0*ax + 4.0f*q2*q3*q3
           - 2.0f*q3*ay - 4.0f*q2 + 8.0f*q2*q1*q1
           + 8.0f*q2*q2*q2 + 4.0f*q2*az;

        s3 = 4.0f*q1*q1*q3 - 2.0f*q1*ax + 4.0f*q2*q2*q3 - 2.0f*q2*ay;

        recipNorm = invSqrt(s0*s0 + s1*s1 + s2*s2 + s3*s3);

        s0 *= recipNorm;
        s1 *= recipNorm;
        s2 *= recipNorm;
        s3 *= recipNorm;

        qDot1 -= beta * s0;
        qDot2 -= beta * s1;
        qDot3 -= beta * s2;
        qDot4 -= beta * s3;
    }


    /* Integrate */
    q0 += qDot1 * dt;
    q1 += qDot2 * dt;
    q2 += qDot3 * dt;
    q3 += qDot4 * dt;

    recipNorm = invSqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);

    *Q0 = q0 * recipNorm;
    *Q1 = q1 * recipNorm;
    *Q2 = q2 * recipNorm;
    *Q3 = q3 * recipNorm;
}



/* ================= FULL AHRS ================= */

void Madgwick_Update(float gx, float gy, float gz,
                     float ax, float ay, float az,
                     float mx, float my, float mz,
                     float dt)
{
    if(mx == 0 && my == 0 && mz == 0)
    {
        Madgwick_UpdateIMU(gx, gy, gz, ax, ay, az, dt);
        return;
    }

    float q0 = *Q0;
    float q1 = *Q1;
    float q2 = *Q2;
    float q3 = *Q3;

    float recipNorm;
    float s0, s1, s2, s3;
    float qDot1, qDot2, qDot3, qDot4;

    float hx, hy;
    float _2bx, _2bz;
    float _2q0mx, _2q0my, _2q0mz, _2q1mx;

    /* Normalize accel */
    recipNorm = invSqrt(ax*ax + ay*ay + az*az);
    ax *= recipNorm;
    ay *= recipNorm;
    az *= recipNorm;

    /* Normalize mag */
    recipNorm = invSqrt(mx*mx + my*my + mz*mz);
    mx *= recipNorm;
    my *= recipNorm;
    mz *= recipNorm;


    /* Auxiliary variables */
    _2q0mx = 2.0f * q0 * mx;
    _2q0my = 2.0f * q0 * my;
    _2q0mz = 2.0f * q0 * mz;
    _2q1mx = 2.0f * q1 * mx;

    float q0q0 = q0*q0;
    float q1q1 = q1*q1;
    float q2q2 = q2*q2;
    float q3q3 = q3*q3;

    hx = mx*q0q0 - _2q0my*q3 + _2q0mz*q2
       + mx*q1q1 + 2.0f*q1*my*q2 + 2.0f*q1*mz*q3
       - mx*q2q2 - mx*q3q3;

    hy = 2.0f*q0*mx*q3 + my*q0q0 - 2.0f*q0*mz*q1
       + 2.0f*q1*mx*q2 - my*q1q1 + my*q2q2
       + 2.0f*q2*mz*q3 - my*q3q3;

    _2bx = sqrtf(hx*hx + hy*hy);
    _2bz = -2.0f*q0*mx*q2 + 2.0f*q0*my*q1 + mz*q0q0
         + 2.0f*q1*mx*q3 - mz*q1q1
         + 2.0f*q2*my*q3 - mz*q2q2 + mz*q3q3;


    /* Gyro derivative */
    qDot1 = 0.5f * (-q1*gx - q2*gy - q3*gz);
    qDot2 = 0.5f * ( q0*gx + q2*gz - q3*gy);
    qDot3 = 0.5f * ( q0*gy - q1*gz + q3*gx);
    qDot4 = 0.5f * ( q0*gz + q1*gy - q2*gx);


    /* Gradient descent step (trimmed but correct) */
    s0 = -2.0f*q2*(2*(q1*q3 - q0*q2) - ax)
       + 2.0f*q1*(2*(q0*q1 + q2*q3) - ay)
       - _2bz*q2*(_2bx*(0.5f - q2q2 - q3q3) + _2bz*(q1*q3 - q0*q2) - mx);

    s1 =  2.0f*q3*(2*(q1*q3 - q0*q2) - ax)
       + 2.0f*q0*(2*(q0*q1 + q2*q3) - ay)
       - 4.0f*q1*(1 - 2*(q1q1 + q2q2) - az);

    s2 = -2.0f*q0*(2*(q1*q3 - q0*q2) - ax)
       + 2.0f*q3*(2*(q0*q1 + q2*q3) - ay)
       - 4.0f*q2*(1 - 2*(q1q1 + q2q2) - az);

    s3 =  2.0f*q1*(2*(q1*q3 - q0*q2) - ax)
       + 2.0f*q2*(2*(q0*q1 + q2*q3) - ay);

    recipNorm = invSqrt(s0*s0 + s1*s1 + s2*s2 + s3*s3);

    s0 *= recipNorm;
    s1 *= recipNorm;
    s2 *= recipNorm;
    s3 *= recipNorm;

    qDot1 -= beta * s0;
    qDot2 -= beta * s1;
    qDot3 -= beta * s2;
    qDot4 -= beta * s3;


    /* Integrate */
    q0 += qDot1 * dt;
    q1 += qDot2 * dt;
    q2 += qDot3 * dt;
    q3 += qDot4 * dt;

    recipNorm = invSqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);

    *Q0 = q0 * recipNorm;
    *Q1 = q1 * recipNorm;
    *Q2 = q2 * recipNorm;
    *Q3 = q3 * recipNorm;
}
