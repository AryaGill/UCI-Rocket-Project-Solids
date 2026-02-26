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
#include "stm32l4xx_hal.h"

/* Private includes ----------------------------------------------------------*/
/* USER CODE BEGIN Includes */

/* USER CODE END Includes */

/* Exported types ------------------------------------------------------------*/
/* USER CODE BEGIN ET */

/* USER CODE END ET */

/* Exported constants --------------------------------------------------------*/
/* USER CODE BEGIN EC */

/* USER CODE END EC */

/* Exported macro ------------------------------------------------------------*/
/* USER CODE BEGIN EM */

/* USER CODE END EM */

void HAL_TIM_MspPostInit(TIM_HandleTypeDef *htim);

/* Exported functions prototypes ---------------------------------------------*/
void Error_Handler(void);

/* USER CODE BEGIN EFP */

/* USER CODE END EFP */

/* Private defines -----------------------------------------------------------*/
#define B1_Pin GPIO_PIN_13
#define B1_GPIO_Port GPIOC
#define RF_CS_Pin GPIO_PIN_0
#define RF_CS_GPIO_Port GPIOC
#define RF_RST_Pin GPIO_PIN_1
#define RF_RST_GPIO_Port GPIOC
#define USART_TX_Pin GPIO_PIN_2
#define USART_TX_GPIO_Port GPIOA
#define USART_RX_Pin GPIO_PIN_3
#define USART_RX_GPIO_Port GPIOA
#define LD2_Pin GPIO_PIN_5
#define LD2_GPIO_Port GPIOA
#define Airbrakes_PWM_1_Pin GPIO_PIN_6
#define Airbrakes_PWM_1_GPIO_Port GPIOA
#define Airbrakes_PWM_2_Pin GPIO_PIN_7
#define Airbrakes_PWM_2_GPIO_Port GPIOA
#define RF_EN_Pin GPIO_PIN_0
#define RF_EN_GPIO_Port GPIOB
#define TMS_Pin GPIO_PIN_13
#define TMS_GPIO_Port GPIOA
#define TCK_Pin GPIO_PIN_14
#define TCK_GPIO_Port GPIOA
#define SWO_Pin GPIO_PIN_3
#define SWO_GPIO_Port GPIOB

/* USER CODE BEGIN Private defines */

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
	float icm_accel_r;
	float icm_accel_p;
	float icm_accel_y;
	float icm_gyro_r;
	float icm_gyro_p;
	float icm_gyro_y;
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
} Telemetry_t;

/* USER CODE END Private defines */

#ifdef __cplusplus
}
#endif

#endif /* __MAIN_H */
