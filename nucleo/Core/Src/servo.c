#include "servo.h"

// Servo timer
extern TIM_HandleTypeDef htim3;

void init_airbrakes_servo(){
	HAL_TIM_PWM_Start(&htim3, AIRBRAKES_SERVO_1_CHANNEL);
	HAL_TIM_PWM_Start(&htim3, AIRBRAKES_SERVO_2_CHANNEL);
}

void set_airbrakes_servo_angle(uint8_t angle)
{
    if (angle > 180) angle = 180;

    uint32_t pulse =
        SERVO_MIN_US +
        ((SERVO_MAX_US - SERVO_MIN_US) * angle) / 180;

    __HAL_TIM_SET_COMPARE(&htim3, AIRBRAKES_SERVO_1_CHANNEL, pulse);
    __HAL_TIM_SET_COMPARE(&htim3, AIRBRAKES_SERVO_2_CHANNEL, pulse);
}
