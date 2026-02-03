#ifndef MADGWICK_FILTER_H
#define MADGWICK_FILTER_H

#ifdef __cplusplus
extern "C" {
#endif

/* Initialization */
void Madgwick_Init(float *q0, float *q1, float *q2, float *q3, float ax, float ay, float az, float b);


/* Full AHRS update (gyro + accel + mag) */
void Madgwick_Update(float gx, float gy, float gz,
                     float ax, float ay, float az,
                     float mx, float my, float mz,
                     float dt);


/* IMU-only update fallback */
void Madgwick_UpdateIMU(float gx, float gy, float gz,
                        float ax, float ay, float az,
                        float dt);

#ifdef __cplusplus
}
#endif

#endif
