#pragma once

#include "main.h"

#define NUM_DEPLOYMENT_LEVELS 128
#define NUM_RECORDED_DEPLOYMENT_LEVELS 11
#define NUM_RECORDED_MACH_NUMS 8

#define SERVO_MIN_US 1000
#define SERVO_MAX_US 2000

#define AIRBRAKES_SERVO_1_CHANNEL TIM_CHANNEL_1
#define AIRBRAKES_SERVO_2_CHANNEL TIM_CHANNEL_2

float get_drag_coefficient(const int deployment_level, const float mach_number);
float get_mach_number(const float velocity, const float temp);
float predict_apogee(Telemetry_t *telemetry, const int deployment_level);
void set_optimal_deployment(FlightState_t flight_state, Telemetry_t *telemetry);
void init_airbrakes_servo();
void set_airbrakes_servo_angle(uint8_t angle);
