#pragma once

#include "main.h"
#include "sd_card.h"

#define FLIGHT_DATA_FILE "flight_data.csv"

void init_data_file();
void get_rf_msg(FlightState_t flight_state, Telemetry_t *t, char* msg, size_t msg_size);
void log_data(FlightState_t flight_state, Telemetry_t *t);
FRESULT write_headers(void);
void state_to_string_num(FlightState_t state, char* str);
void state_to_string_name(FlightState_t state, char* str);
void save_data_file();
void write_datafile_message(char* msg);
uint32_t telemetry_log_period(FlightState_t state);
