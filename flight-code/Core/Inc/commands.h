#pragma once

#include "fsm.h"

#define RX_BUF_SIZE 64

void handle_rf_command(char *cmd, FlightState_t *flight_state, Telemetry_t *telemetry);
