#pragma once

#include "main.h"

#define SERVO_MIN_US 1000
#define SERVO_MAX_US 2000

#define AIRBRAKES_SERVO_1_CHANNEL TIM_CHANNEL_1
#define AIRBRAKES_SERVO_2_CHANNEL TIM_CHANNEL_2

void init_airbrakes_servo();
void set_airbrakes_servo_angle(uint8_t angle);
