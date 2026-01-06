#pragma once

#include "telemetry.h"

void init_sensors();
void read_sensors(Telemetry *telemetry);
void read_bmp(Telemetry *telemetry);
