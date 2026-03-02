#pragma once

#include <string.h>
#include "main.h"

// Change before flight - check airbrakes.h as well
#define MAIN_DEPLOY_MAX_ALT 229 // m
#define MAIN_DEPLOY_MIN_ALT 77 // m
// Leilani rocket
#define MOTOR_BURN_TIME 1500 // ms
#define POWER_RESET_MIN_ALT_CHANGE 15
// Night Fury
//#define MOTOR_BURN_TIME 4000 // ms
//#define POWER_RESET_MIN_ALT_CHANGE 40
// Light Fury
//#define MOTOR_BURN_TIME 4000 // ms
//#define POWER_RESET_MIN_ALT_CHANGE 40

#define LAUNCH_THRESHOLD 10
#define APOGEE_THRESHOLD -0.5
#define LANDED_THRESHOLD -0.2
#define ALT_DIF_BUF_SIZE 10

// Liftoff detection constants
#define LAUNCH_ACCEL_THRESHOLD 40
#define RAIL_DELAY_TIME 250
#define LAUNCH_EVAL_PERIOD_TIME 250

#define STATE_FILE "flight_state.csv"
#define MIN_RESET_ALT 150 // m

float get_avg_alt_dif();
void update_alt_dif_buf(float new_alt_dif);
void set_flight_state(FlightState_t new_state, FlightState_t *flight_state, Telemetry_t *telemetry);
void init_flight_state(FlightState_t *flight_state, Telemetry_t *telemetry);
void update_flight_state(FlightState_t *flight_state, Telemetry_t *telemetry);
