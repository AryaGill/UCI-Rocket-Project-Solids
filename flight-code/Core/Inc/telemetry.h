#pragma once

#include "main.h"
#include "sd_card.h"

#define FLIGHT_DATA_FILE "flight_data.csv"

void init_data_file();
void log_data(char* state_str, Telemetry_t *telemetry);
FRESULT write_headers(void);
