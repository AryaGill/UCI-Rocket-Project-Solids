#pragma once

#include "fsm.h"
#include "telemetry.h"

typedef struct {
	float pressure;
	float altitude;
	float startAlt;
	float temperature;
	float angle_of_attack;
	float velocity_r;
	float velocity_p;
	float velocity_y;
	float accel_r;
	float accel_p;
	float accel_y;
	float gyro_r;
	float gyro_p;
	float gyro_y;
	float predicted_apogee;
	float airbrake_deployment;
	float mag_r;
	float mag_p;
	float mag_y;
} Telemetry_t;

void init_dataFile();
void log_data(FlightState_t *flight_state, Telemetry_t *telemetry);
