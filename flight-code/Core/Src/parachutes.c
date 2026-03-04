#include "parachutes.h"

extern ADC_HandleTypeDef hadc1;

static ADC_ChannelConfTypeDef baseConfig = {
    .Rank = ADC_REGULAR_RANK_1,
    .SamplingTime = ADC_SAMPLETIME_64CYCLES_5,
    .SingleDiff = ADC_SINGLE_ENDED,
    .OffsetNumber = ADC_OFFSET_NONE,
    .Offset = 0,
	.OffsetSignedSaturation = DISABLE
};

uint32_t read_ematch_voltage(ADC_HandleTypeDef *hadc, uint32_t channel){
	ADC_ChannelConfTypeDef sConfig = {0};
	  sConfig.Rank = ADC_REGULAR_RANK_1;
	  sConfig.SamplingTime = ADC_SAMPLETIME_64CYCLES_5;
	  sConfig.SingleDiff = ADC_SINGLE_ENDED;
	  sConfig.OffsetNumber = ADC_OFFSET_NONE;
	  sConfig.Offset = 0;
	  sConfig.OffsetSignedSaturation = DISABLE;

//	ADC_ChannelConfTypeDef sConfig = baseConfig;
	sConfig.Channel = channel;

	if (HAL_ADC_ConfigChannel(hadc, &sConfig) != HAL_OK) return 0xFFFFFFFF;
	if (HAL_ADC_Start(hadc) != HAL_OK) return 0xFFFFFFFF;
	if (HAL_ADC_PollForConversion(hadc, 5) != HAL_OK) return 0xFFFFFFFF;

	uint32_t adc = HAL_ADC_GetValue(hadc);

	HAL_ADC_Stop(hadc);

	return adc;
}

void read_ematch_connections(Telemetry_t *telemetry){
	telemetry->main_p_ematch_voltage = read_ematch_voltage(&hadc1, MAIN_P_ADC_CHANNEL);
	telemetry->main_s_ematch_voltage = read_ematch_voltage(&hadc1, MAIN_S_ADC_CHANNEL);
	telemetry->drogue_p_ematch_voltage = read_ematch_voltage(&hadc1, DROGUE_P_ADC_CHANNEL);
	telemetry->drogue_s_ematch_voltage = read_ematch_voltage(&hadc1, DROGUE_S_ADC_CHANNEL);
}

void drogue_primary_on(){
	HAL_GPIO_WritePin(Drogue_Parachute_1_GPIO_Port, Drogue_Parachute_1_Pin, GPIO_PIN_SET);
}

void drogue_primary_off(){
	HAL_GPIO_WritePin(Drogue_Parachute_1_GPIO_Port, Drogue_Parachute_1_Pin, GPIO_PIN_RESET);
}

void drogue_secondary_on(){
	HAL_GPIO_WritePin(Drogue_Parachute_2_GPIO_Port, Drogue_Parachute_2_Pin, GPIO_PIN_SET);
}

void drogue_secondary_off(){
	HAL_GPIO_WritePin(Drogue_Parachute_2_GPIO_Port, Drogue_Parachute_2_Pin, GPIO_PIN_RESET);
}

void main_primary_on(){
	HAL_GPIO_WritePin(Main_Parachute_1_GPIO_Port, Main_Parachute_1_Pin, GPIO_PIN_SET);
}

void main_primary_off(){
	HAL_GPIO_WritePin(Main_Parachute_1_GPIO_Port, Main_Parachute_1_Pin, GPIO_PIN_RESET);
}

void main_secondary_on(){
	HAL_GPIO_WritePin(Main_Parachute_2_GPIO_Port, Main_Parachute_2_Pin, GPIO_PIN_SET);
}

void main_secondary_off(){
	HAL_GPIO_WritePin(Main_Parachute_2_GPIO_Port, Main_Parachute_2_Pin, GPIO_PIN_RESET);
}
