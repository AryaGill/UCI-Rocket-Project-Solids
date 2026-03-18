#ifndef PARACHUTE_H
#define PARACHUTE_H

#include <stdbool.h>
#include <stdint.h>
#include "main.h"

// Velocity bands (ft/s) for baro-based validation
#define PARACHUTE_DROGUE_VEL_MIN      (-22.08f)   // faster (more negative)
#define PARACHUTE_DROGUE_VEL_MAX      (-19.98f)   // slower (less negative)

#define PARACHUTE_MAIN_VEL_TARGET    (-8.97f)


// Sample counts for confirmation (tune vs. sample rate)
#define PARACHUTE_DROGUE_SAMPLES_REQUIRED  10
#define PARACHUTE_MAIN_SAMPLES_REQUIRED    10


/**
 * @brief Update baro-based validation flags for drogue and main parachutes.
 *
 * Expects:
 *  - Telemetry_t has vertical velocity in ft/s (negative = descending),
 *    e.g. `float baro_vel_ftps;`.
 *  - FlightState_t contains at least:
 *        bool apogee_detected;
 *        bool drogue_fired;
 *        bool main_fired;
 *        bool drogue_validated_baro;
 *        bool main_validated_baro;
 *
 * Call once per main loop iteration after filters and flight_state update.
 */
void parachute_update_recovery_validation(FlightState_t *state, Telemetry_t *t);

#endif // PARACHUTE_H
