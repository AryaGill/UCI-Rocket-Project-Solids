//#pragma once
//#include "main.h"

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
#define TAU_BARO_VEL   0.5f // seconds

#define BARO_TRUST_TIME 4.0f

//Initialize and update complementary filter
void complementary_filter_init(Telemetry_t *telemetry);
void complementary_filter(Telemetry_t *telemetry, FlightState_t *flight_state);

#ifdef __cplusplus
}
#endif

#endif /* COMPLEMENTARY_FILTER_H */
