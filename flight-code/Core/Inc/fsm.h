#pragma once

#include <string.h>

#define LAUNCH_THRESHOLD 10
#define APOGEE_THRESHOLD -0.5
#define LANDED_THRESHOLD -0.2
#define ALT_DIF_BUF_SIZE 10

typedef enum {
  LAUNCH_PAD,
  MOTOR_BURN,
  GLIDING_ASCENT,
  DROGUE_PRIMARY_DEPLOYING,
  DROGUE_PRIMARY_DEPLOYED,
  DROGUE_SECONDARY_DEPLOYING,
  DROGUE_SECONDARY_DEPLOYED,
  MAIN_PRIMARY_DEPLOYING,
  MAIN_PRIMARY_DEPLOYED,
  MAIN_SECONDARY_DEPLOYING,
  MAIN_SECONDARY_DEPLOYED,
  LANDED
} FlightState_t;

float get_avg_alt_dif();
void update_alt_dif_buf(float new_alt_dif);
void set_flight_state(FlightState_t new_state, FlightState_t *flight_state);
void initialize_flight_state(FlightState_t *flight_state, Telemetry_t *telemetry);
void update_flight_state(FlightState_t *flight_state, Telemetry_t *telemetry);
String state_to_string(FlightState state);
