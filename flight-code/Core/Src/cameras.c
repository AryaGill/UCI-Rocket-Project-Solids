#include "cameras.h"

extern ADC_HandleTypeDef hadc1;

GPIO_TypeDef *CAM_PORTS[2] = {Camera_1_GPIO_Port, Camera_2_GPIO_Port};
uint16_t CAM_PINS[2] = {Camera_1_Pin, Camera_2_Pin};
ADC_HandleTypeDef *CAM_ADCs[2] = {&hadc1, &hadc1};
uint32_t CAM_ADC_CHANNELS[2] = {ADC_CHANNEL_9, ADC_CHANNEL_5};

void turn_camera_on(int cam_num){
	HAL_GPIO_WritePin(CAM_PORTS[cam_num], CAM_PINS[cam_num], GPIO_PIN_SET);
//	HAL_Delay(5);
//
//	ADC_ChannelConfTypeDef sConfig;
//	for(int i = 0; i < 5; ++i){
//		// Read ADC
//		sConfig.Channel = CAM_ADC_CHANNELS[cam_num];
//		HAL_ADC_ConfigChannel(CAM_ADCs[cam_num], &sConfig);
//
//		HAL_ADC_Start(CAM_ADCs[cam_num]);
//		HAL_ADC_PollForConversion(CAM_ADCs[cam_num], HAL_MAX_DELAY);
//
//		uint32_t adc = HAL_ADC_GetValue(CAM_ADCs[cam_num]);
//
//		HAL_ADC_Stop(CAM_ADCs[cam_num]);
//
//		if (adc < 1000){
//			// Camera is on
//			break;
//		}
//		else{
//			// Camera is off. Turn off then back on.
//			HAL_GPIO_WritePin(CAM_PORTS[cam_num], CAM_PINS[cam_num], GPIO_PIN_RESET);
//			HAL_Delay(10);
//			HAL_GPIO_WritePin(CAM_PORTS[cam_num], CAM_PINS[cam_num], GPIO_PIN_SET);
//			HAL_Delay(10);
//		}
//	}
}

void turn_camera_off(int cam_num){
	HAL_GPIO_WritePin(CAM_PORTS[cam_num], CAM_PINS[cam_num], GPIO_PIN_RESET);
}

uint32_t read_adc(uint32_t channel, ADC_HandleTypeDef *hadc){
	ADC_ChannelConfTypeDef sConfig = {0};
	sConfig.Rank = ADC_REGULAR_RANK_1;
	sConfig.SamplingTime = ADC_SAMPLETIME_64CYCLES_5;
	sConfig.SingleDiff = ADC_SINGLE_ENDED;
	sConfig.OffsetNumber = ADC_OFFSET_NONE;
	sConfig.Offset = 0;
	sConfig.OffsetSignedSaturation = DISABLE;

	sConfig.Channel = channel;
	if (HAL_ADC_ConfigChannel(hadc, &sConfig) != HAL_OK) return 0;
	if (HAL_ADC_Start(hadc) != HAL_OK) return 0;
	if (HAL_ADC_PollForConversion(hadc, 5) != HAL_OK) return 0;

	uint32_t adc = HAL_ADC_GetValue(hadc);

	HAL_ADC_Stop(hadc);

	return adc;
}

uint32_t cam1_adc;
uint32_t cam2_adc;
void read_camera_adcs(Telemetry_t *telemetry){
	// Read ADC 1
	cam1_adc = read_adc(CAM_ADC_CHANNELS[0], CAM_ADCs[0]);
	telemetry->cam1_on = (cam1_adc < 1000) ? 1 : 0;

	// Read ADC 2
	cam2_adc = read_adc(CAM_ADC_CHANNELS[1], CAM_ADCs[1]);
	telemetry->cam2_on = (cam2_adc < 1000) ? 1 : 0;
}
