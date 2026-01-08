#pragma once

#include "fsm.h"
#include "telemetry.h"

#define NUM_DEPLOYMENT_LEVELS 1024
#define NUM_RECORDED_DEPLOYMENT_LEVELS 11
#define NUM_RECORDED_MACH_NUMS 8

// Deployment levels should be evenly spread between least and most deployment (inclusive)
// Mach numbers should be evenly spread between 0 and 0.7 (inclusive)
float air_brakes_drag_coefficient[NUM_RECORDED_DEPLOYMENT_LEVELS][NUM_RECORDED_MACH_NUMS] = {
    {0.400, 0.400, 0.400, 0.400, 0.400, 0.400, 0.400, 0.400}, // deploy 0.0
    {0.430, 0.460, 0.490, 0.520, 0.550, 0.580, 0.610, 0.640}, // deploy 0.1
    {0.460, 0.490, 0.520, 0.550, 0.580, 0.610, 0.640, 0.670}, // deploy 0.2
    {0.490, 0.520, 0.550, 0.580, 0.610, 0.640, 0.670, 0.700}, // deploy 0.3
    {0.520, 0.550, 0.580, 0.610, 0.640, 0.670, 0.700, 0.730}, // deploy 0.4
    {0.550, 0.580, 0.610, 0.640, 0.670, 0.700, 0.730, 0.760}, // deploy 0.5
    {0.580, 0.610, 0.640, 0.670, 0.700, 0.730, 0.760, 0.790}, // deploy 0.6
    {0.610, 0.640, 0.670, 0.700, 0.730, 0.760, 0.790, 0.820}, // deploy 0.7
    {0.640, 0.670, 0.700, 0.730, 0.760, 0.790, 0.820, 0.850}, // deploy 0.8
    {0.670, 0.700, 0.730, 0.760, 0.790, 0.820, 0.850, 0.880}, // deploy 0.9
    {0.700, 0.730, 0.760, 0.790, 0.820, 0.850, 0.880, 0.910}  // deploy 1.0
};

float get_drag_coefficient(const int deployment_level, const float mach_number);
float get_mach_number(const float velocity, const float temp);
float predict_apogee(Telemetry_t *telemetry, const int deployment_level);
void set_optimal_deployment(FlightState_t flight_state, Telemetry_t *telemetry);
