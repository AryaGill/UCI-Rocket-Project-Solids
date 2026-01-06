#pragma once

#include "fsm.h"
#include "telemetry.h"

typedef struct {
	float altitude;
	float startAlt;
	float temperature;
} Telemetry_t;

void init_dataFile();
void log_data(FlightState_t *flight_state, Telemetry_t *telemetry);
