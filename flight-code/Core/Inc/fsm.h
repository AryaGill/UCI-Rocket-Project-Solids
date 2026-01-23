#pragma once

#include <string.h>
#include "main.h"

#define LAUNCH_THRESHOLD 10
#define APOGEE_THRESHOLD -0.5
#define LANDED_THRESHOLD -0.2
#define ALT_DIF_BUF_SIZE 10

#define MOTOR_BURN_TIME 1500

#define MAIN_DEPLOY_MAX_ALT 229
#define MAIN_DEPLOY_MIN_ALT 77

float get_avg_alt_dif();
void update_alt_dif_buf(float new_alt_dif);
void set_flight_state(FlightState_t new_state, FlightState_t *flight_state);
void init_flight_state(FlightState_t *flight_state, Telemetry_t *telemetry);
void update_flight_state(FlightState_t *flight_state, Telemetry_t *telemetry);
void state_to_string(FlightState_t state, char* str);
