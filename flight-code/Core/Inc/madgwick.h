//#pragma once
//#include "main.h"
//
///* Initialization */
//void Madgwick_Init(Telemetry_t* telemetry, float b);
//
///* Full AHRS update (gyro + accel + mag) */
//void Madgwick_Update(Telemetry_t* telemetry);
//
///* IMU-only update fallback */
//void Madgwick_UpdateIMU(Telemetry_t* telemetry);
//
//

#ifndef MADGWICK_H
#define MADGWICK_H

#include "telemetry.h"
#include <stdint.h>

void Madgwick_Init(Telemetry_t* telemetry, float beta);
void Madgwick_Update(Telemetry_t* telemetry);
void Madgwick_GetEuler(Telemetry_t* telemetry);

#endif // MADGWICK_H
