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
	float predicted_apogee;
	float airbrake_deployment;
} Telemetry_t;

void init_dataFile();
void log_data(FlightState_t *flight_state, Telemetry_t *telemetry);
