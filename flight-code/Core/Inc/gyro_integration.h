#pragma once
#include "main.h"

typedef struct{
	float q0;
	float q1;
	float q2;
	float q3;
} Quaternions_t;

void Gyro_Integrate(Telemetry_t* telemetry, Quaternions_t* quats);
void Quaternion_From_Accel(Telemetry_t* telemetry, Quaternions_t* quats);
void integrate_gyro(FlightState_t flight_state, Telemetry_t* telemetry);
