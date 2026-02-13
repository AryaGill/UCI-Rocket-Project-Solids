#pragma once
#include "main.h"

#define EMATCH_CONNECTED_THRESHOLD 100

void turn_camera_on(int cam_num);
void turn_camera_off(int cam_num);
void read_camera_adcs(Telemetry_t *telemetry);
