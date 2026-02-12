#pragma once
#include "main.h"

extern ADC_HandleTypeDef hadc1;

void turn_camera_on(int cam_num);
void turn_camera_off(int cam_num);
void read_camera_adcs(Telemetry_t *telemetry);
