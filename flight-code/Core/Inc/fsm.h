#pragma once

#include <string.h>
#include "main.h"

// Change before flight - check airbrakes.h as well
#define MAIN_DEPLOY_MAX_ALT 304 // m (1000 ft)
#define MAIN_DEPLOY_MIN_ALT 50 // m
#define DROGUE_DEPLOY_MIN_ALT 50 // m
// Leilani rocket
#define MOTOR_BURN_TIME 2800 // ms
// Night Fury
// #define MOTOR_BURN_TIME 5000 // ms
// Bright Fury
// #define MOTOR_BURN_TIME 6000` // ms

#define POWER_RESET_MIN_ALT_CHANGE 10

#define APOGEE_VELO_THRESHOLD -1.0
#define LANDED_VELO_THRESHOLD -0.2

// Liftoff detection constants
#define LAUNCH_ACCEL_THRESHOLD 40
#define RAIL_DELAY_TIME 250
#define LAUNCH_EVAL_PERIOD_TIME 250// TODO: set back to 250

#define STATE_FILE "flight_state.csv"
#define MIN_RESET_ALT 50 // m

//float get_avg_alt_dif();
//void update_alt_dif_buf(float new_alt_dif);
void set_flight_state(FlightState_t new_state, FlightState_t *flight_state, Telemetry_t *telemetry);
void init_flight_state(FlightState_t *flight_state, Telemetry_t *telemetry);
void update_flight_state(FlightState_t *flight_state, Telemetry_t *telemetry);
