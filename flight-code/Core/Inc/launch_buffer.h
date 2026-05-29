#ifndef LAUNCH_BUFFER_H
#define LAUNCH_BUFFER_H

#include <stdint.h>
#include "telemetry.h"
#include "fsm.h"

//number of entries stored in buffer
#define LAUNCH_BUFFER_SIZE 200

//max num characters within each entry
#define LAUNCH_BUFFER_LINE_SIZE 1000

//implement circular buffer
//head: stores the oldest entry
//full: reset head if full or continue filling buffer till full
//flushed: if all data is flushed already stop all buffer operations
typedef struct {
    char buffer[LAUNCH_BUFFER_SIZE][LAUNCH_BUFFER_LINE_SIZE];
    uint16_t head;
    uint8_t full;
    uint8_t flushed;
} launch_buffer_t;

//initialize buffer, add entry to buffer, flush all of buffer into SD card
void launch_buffer_init(launch_buffer_t *lb);
void launch_buffer_add(launch_buffer_t *lb, FlightState_t flight_state, Telemetry_t *t);
void launch_buffer_flush(launch_buffer_t *lb);

#endif
