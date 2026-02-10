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
	float pressure;
	float altitude;
	float startAlt;
	float temperature;
	float angle_of_attack;
	float velocity_r;
	float velocity_p;
	float velocity_y;
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
	float airbrake_deployment;
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
} Telemetry_t;

/* USER CODE END ET */

/* Exported constants --------------------------------------------------------*/
/* USER CODE BEGIN EC */

/* USER CODE END EC */

/* Exported macro ------------------------------------------------------------*/
/* USER CODE BEGIN EM */

/* USER CODE END EM */

/* Exported functions prototypes ---------------------------------------------*/
void Error_Handler(void);

/* USER CODE BEGIN EFP */

uint32_t micros(void);

/* USER CODE END EFP */

/* Private defines -----------------------------------------------------------*/
#define LED_Pin GPIO_PIN_0
#define LED_GPIO_Port GPIOC
#define Flash_CS_Pin GPIO_PIN_3
#define Flash_CS_GPIO_Port GPIOA
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
#define Drogue_Parachute_2_Pin GPIO_PIN_10
#define Drogue_Parachute_2_GPIO_Port GPIOD
#define Drogue_Parachute_1_Pin GPIO_PIN_11
#define Drogue_Parachute_1_GPIO_Port GPIOD
#define Mag_SDIO_Pin GPIO_PIN_6
#define Mag_SDIO_GPIO_Port GPIOC
#define Mag_CS_Pin GPIO_PIN_7
#define Mag_CS_GPIO_Port GPIOC
#define Main_Parachute_2_Pin GPIO_PIN_9
#define Main_Parachute_2_GPIO_Port GPIOA
#define Main_Parachute_1_Pin GPIO_PIN_10
#define Main_Parachute_1_GPIO_Port GPIOA
#define Buzzer_Pin GPIO_PIN_0
#define Buzzer_GPIO_Port GPIOD

/* USER CODE BEGIN Private defines */

#define CHARGE_DELAY 500
#define BACKUP_DELAY 500

/* USER CODE END Private defines */

#ifdef __cplusplus
}
#endif

#endif /* __MAIN_H */
