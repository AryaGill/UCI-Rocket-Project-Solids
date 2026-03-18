/* USER CODE BEGIN Header */
/**
  ******************************************************************************
  * @file           : main.h
  * @brief          : Header for main.c file.
  *                   This file contains the common defines of the application.
  ******************************************************************************
  * @attention
  *
  * Copyright (c) 2026 STMicroelectronics.
  * All rights reserved.
  *
  * This software is licensed under terms that can be found in the LICENSE file
  * in the root directory of this software component.
  * If no LICENSE file comes with this software, it is provided AS-IS.
  *
  ******************************************************************************
  */
/* USER CODE END Header */

/* Define to prevent recursive inclusion -------------------------------------*/
#ifndef __MAIN_H
#define __MAIN_H

#ifdef __cplusplus
extern "C" {
#endif

/* Includes ------------------------------------------------------------------*/
#include "stm32h7xx_hal.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */

/* USER CODE END Includes */

/* Exported types ------------------------------------------------------------*/
/* USER CODE BEGIN ET */

typedef enum {
	DISARMED,
	LAUNCH_PAD,
	MOTOR_BURN,
	GLIDING_ASCENT,
	DROGUE_PRIMARY_DEPLOYING,
	DROGUE_PRIMARY_DEPLOYED,
	DROGUE_SECONDARY_DEPLOYING,
	DROGUE_SECONDARY_DEPLOYED,
	MAIN_PRIMARY_DEPLOYING,
	MAIN_PRIMARY_DEPLOYED,
	MAIN_SECONDARY_DEPLOYING,
	MAIN_SECONDARY_DEPLOYED,
	LANDED
} FlightState_t;

typedef struct {
	uint32_t time;
	float pressure;
	float altitude;
	float startAlt;
	float temperature;
//	float angle_of_attack;
	float velocity_world_x;
	float velocity_world_y;
	float velocity_world_z;
//	float velocity_r;
//	float velocity_p;
//	float velocity_y;
	float lsm_accel_r;
	float lsm_accel_p;
	float lsm_accel_y;
	float lsm_gyro_r;
	float lsm_gyro_p;
	float lsm_gyro_y;
	float adxl_accel_r;
	float adxl_accel_p;
	float adxl_accel_y;
	float predicted_apogee;
	uint8_t airbrake_deployment;
	float mag_r;
	float mag_p;
	float mag_y;
	float q0;
	float q1;
	float q2;
	float q3;
	float accel_world_x;
	float accel_world_y;
	float accel_world_z;
	float alt_fused;
	uint8_t cam1_on;
	uint8_t cam2_on;
	uint32_t main_p_ematch_voltage;
	uint32_t main_s_ematch_voltage;
	uint32_t drogue_p_ematch_voltage;
	uint32_t drogue_s_ematch_voltage;
	float roll;
	float pitch;
	float yaw;
	uint32_t t_burnout;
	uint32_t t_apogee;
	uint32_t t_drogue;
	uint32_t t_main;
	uint32_t t_land;
	float baro_vz;
	int drogue_validated_baro;
	int main_validated_baro;
} Telemetry_t;

typedef struct{
	float adxl_accel_r_bias;
	float adxl_accel_p_bias;
	float adxl_accel_y_bias;

	float lsm_accel_r_bias;
	float lsm_accel_p_bias;
	float lsm_accel_y_bias;

	float lsm_gyro_r_bias;
	float lsm_gyro_p_bias;
	float lsm_gyro_y_bias;

	float mag_r_bias;
	float mag_p_bias;
	float mag_y_bias;

	float mag_r_scale;
	float mag_p_scale;
	float mag_y_scale;

	uint32_t bias_count;
} Bias_t;
/* USER CODE END ET */

/* Exported constants --------------------------------------------------------*/
/* USER CODE BEGIN EC */

#define RF_TRANSMIT_PERIOD 500 // ms

/* USER CODE END EC */

/* Exported macro ------------------------------------------------------------*/
/* USER CODE BEGIN EM */

/* USER CODE END EM */

void HAL_TIM_MspPostInit(TIM_HandleTypeDef *htim);

/* Exported functions prototypes ---------------------------------------------*/
void Error_Handler(void);

/* USER CODE BEGIN EFP */

uint64_t micros(void);

/* USER CODE END EFP */

/* Private defines -----------------------------------------------------------*/
#define RF_EN_Pin GPIO_PIN_2
#define RF_EN_GPIO_Port GPIOE
#define RF_CS_Pin GPIO_PIN_3
#define RF_CS_GPIO_Port GPIOE
#define RF_RST_Pin GPIO_PIN_4
#define RF_RST_GPIO_Port GPIOE
#define LED_Pin GPIO_PIN_0
#define LED_GPIO_Port GPIOC
#define RX_RF_Pin GPIO_PIN_0
#define RX_RF_GPIO_Port GPIOA
#define RF_TX_Pin GPIO_PIN_1
#define RF_TX_GPIO_Port GPIOA
#define Drogue_ADC1_Pin GPIO_PIN_2
#define Drogue_ADC1_GPIO_Port GPIOA
#define Flash_CS_Pin GPIO_PIN_3
#define Flash_CS_GPIO_Port GPIOA
#define Drogue_ADC2_Pin GPIO_PIN_4
#define Drogue_ADC2_GPIO_Port GPIOA
#define Main_ADC1_Pin GPIO_PIN_4
#define Main_ADC1_GPIO_Port GPIOC
#define Main_ADC2_Pin GPIO_PIN_5
#define Main_ADC2_GPIO_Port GPIOC
#define CAM_ADC1_Pin GPIO_PIN_0
#define CAM_ADC1_GPIO_Port GPIOB
#define CAM_ADC2_Pin GPIO_PIN_1
#define CAM_ADC2_GPIO_Port GPIOB
#define IMU_2_CS_Pin GPIO_PIN_2
#define IMU_2_CS_GPIO_Port GPIOB
#define SD_CS_Pin GPIO_PIN_9
#define SD_CS_GPIO_Port GPIOE
#define SD_CD_Pin GPIO_PIN_10
#define SD_CD_GPIO_Port GPIOE
#define Baro_CS_Pin GPIO_PIN_11
#define Baro_CS_GPIO_Port GPIOE
#define IMU_CS_Pin GPIO_PIN_10
#define IMU_CS_GPIO_Port GPIOB
#define Buzzer_Pin GPIO_PIN_12
#define Buzzer_GPIO_Port GPIOD
#define Camera_1_Pin GPIO_PIN_14
#define Camera_1_GPIO_Port GPIOD
#define Camera_2_Pin GPIO_PIN_15
#define Camera_2_GPIO_Port GPIOD
#define Mag_CS_Pin GPIO_PIN_7
#define Mag_CS_GPIO_Port GPIOC
#define Drogue_Parachute_2_Pin GPIO_PIN_9
#define Drogue_Parachute_2_GPIO_Port GPIOC
#define Drogue_Parachute_1_Pin GPIO_PIN_8
#define Drogue_Parachute_1_GPIO_Port GPIOA
#define Main_Parachute_2_Pin GPIO_PIN_9
#define Main_Parachute_2_GPIO_Port GPIOA
#define Main_Parachute_1_Pin GPIO_PIN_10
#define Main_Parachute_1_GPIO_Port GPIOA
#define Airbrakes_PWM_1_Pin GPIO_PIN_4
#define Airbrakes_PWM_1_GPIO_Port GPIOB
#define Airbrakes_PWM_2_Pin GPIO_PIN_5
#define Airbrakes_PWM_2_GPIO_Port GPIOB

/* USER CODE BEGIN Private defines */

#define CHARGE_DELAY 500
#define BACKUP_DELAY 500

/* USER CODE END Private defines */

#ifdef __cplusplus
}
#endif

#endif /* __MAIN_H */
