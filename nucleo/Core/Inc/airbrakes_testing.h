#pragma once

#include "main.h"

#define NUM_DEPLOYMENT_LEVELS 64
#define NUM_RECORDED_DEPLOYMENT_LEVELS 11
#define NUM_RECORDED_MACH_NUMS 14

#define SERVO_MIN_US 1000
#define SERVO_MAX_US 2000

#define AIRBRAKES_SERVO_1_CHANNEL TIM_CHANNEL_1
#define AIRBRAKES_SERVO_2_CHANNEL TIM_CHANNEL_2

float get_CdA(uint8_t deployment_level, const float mach_number);
float get_mach_number(const float velocity, const float temp);
float predict_apogee(Telemetry_t *telemetry, uint8_t deployment_level);
void set_optimal_deployment(FlightState_t flight_state, Telemetry_t *telemetry);
void init_airbrakes_servo();
void set_airbrakes_servo_angle(uint8_t angle);
void set_airbrakes_deployment_level(Telemetry_t *telemetry, uint8_t deployment);
