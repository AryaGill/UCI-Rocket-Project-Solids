#pragma once
#include "main.h"

extern ADC_HandleTypeDef hadc1;

GPIO_TypeDef *CAM_PORTS[2] = {Camera_1_GPIO_Port, Camera_2_GPIO_Port};
uint16_t CAM_PINS[2] = {Camera_1_Pin, Camera_2_Pin};
ADC_HandleTypeDef *CAM_ADCs[2] = {&hadc1, &hadc1};
uint32_t CAM_ADC_CHANNELS[2] = {ADC_CHANNEL_9, ADC_CHANNEL_5};

void turn_camera_on(int cam_num);
void turn_camera_off(int cam_num);
