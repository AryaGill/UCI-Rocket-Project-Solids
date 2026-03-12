#pragma once

#include "main.h"

// Change before flight - check fsm.h as well
#define MASS 23.77 // kg
#define TARGET_APOGEE_FT 8000
#define TARGET_APOGEE_M TARGET_APOGEE_FT * 0.3048

#define SERVO_ANGLE_NOT_EXTENDED 180.0f
#define SERVO_ANGLE_EXTENDED 137.0f

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
void set_airbrakes_servo_angle(float angle);
void set_airbrakes_deployment_level(Telemetry_t *telemetry, uint8_t deployment);
void set_airbrakes_initial_temp(Telemetry_t *telemetry);
void perform_airbrakes_servo_sequence();
