#pragma once

#include "fsm.h"
#include "telemetry.h"

void set_airbrake_deployment(int deployment);
void optimal_deployment(FlightState_t *flight_state, Telemetry_t *telemetry); // Make sure only works when flight_state == GLIDING_ASCENT && get_mach_number(velocity, Temp) < 0.7 && angle of attack < 30 deg... else 0
