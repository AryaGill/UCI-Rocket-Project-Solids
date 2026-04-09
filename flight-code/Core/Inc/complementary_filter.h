//#pragma once
//#include "main.h"
//
//// Complementary filter weights
//#define ALPHA_VELOCITY 0.99f  // 99% IMU integrated velocity
//#define ALPHA_ALTITUDE 0.95f  // 95% integrated fused velocity
//
//void complementary_filter(Telemetry_t* telemetry);

#ifndef COMPLEMENTARY_FILTER_H
#define COMPLEMENTARY_FILTER_H

#include "telemetry.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Complementary filter gains
 *
 * ALPHA_VELOCITY:
 *   Higher = trust IMU-integrated vertical velocity more
 *   Lower  = trust baro-derived vertical velocity more
 *
 * ALPHA_ALTITUDE:
 *   Higher = trust integrated altitude prediction more
 *   Lower  = trust raw barometer altitude more
 *
 * Both should stay in the range [0.0, 1.0].
 */
#define TAU_VELOCITY   0.2f // seconds
#define TAU_ALTITUDE   0.2f // seconds
#define TAU_BARO_VEL   0.4f // seconds

void complementary_filter_init(Telemetry_t *telemetry);
void complementary_filter(Telemetry_t *telemetry, FlightState_t *flight_state);

#ifdef __cplusplus
}
#endif

#endif /* COMPLEMENTARY_FILTER_H */
