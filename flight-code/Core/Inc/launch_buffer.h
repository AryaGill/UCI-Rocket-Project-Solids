#ifndef LAUNCH_BUFFER_H
#define LAUNCH_BUFFER_H

#include <stdint.h>
#include "telemetry.h"
#include "fsm.h"

#define LAUNCH_BUFFER_SIZE 200
#define LAUNCH_BUFFER_LINE_SIZE 1000

typedef struct {
    char buffer[LAUNCH_BUFFER_SIZE][LAUNCH_BUFFER_LINE_SIZE];
    uint16_t head;
    uint8_t full;
    uint8_t flushed;
} launch_buffer_t;

void launch_buffer_init(launch_buffer_t *lb);
void launch_buffer_add(launch_buffer_t *lb, FlightState_t flight_state, Telemetry_t *t);
void launch_buffer_flush(launch_buffer_t *lb);

#endif
