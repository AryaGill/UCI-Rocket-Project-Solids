#pragma once

#include "main.h"

// LSMDSL (IMU 1)
#define LSM6DSL_CTRL1_XL   0x10
#define LSM6DSL_CTRL2_G    0x11
#define LSM6DSL_CTRL3_C    0x12
#define LSM6DSL_OUTX_L_G   0x22
#define LSM6DSL_WHO_AM_I 0x0F //temp remove later
// ICM45686 (IMU 2)
#define ICM_PWR_MGMT0   0x4E
#define ICM_GYRO_CFG0   0x4F
#define ICM_ACCEL_CFG0  0x50
#define ICM_DATA_START  0x1F

// LPS22HHTR (Baro)
#define LPS22HH_CTRL_REG1  0x10
#define LPS22HH_PRESS_OUT  0x28
#define LPS22HH_TEMP_OUT   0x2B

// IIS2MDCTR (Mag)
#define IIS2M_CTRL_REG1  0x20
#define IIS2M_OUTX_L    0x28

void init_sensors(SPI_HandleTypeDef *hspi);
void read_sensors(Telemetry_t *telemetry);
uint8_t LSM6DSL_WhoAmI(void);
void LSM6DSL_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin);
void LSM6DSL_Read(Telemetry_t *telemetry);
void ICM45686_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin);
void ICM45686_Read(Telemetry_t *telemetry);
void LPS22HH_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin);
void LPS22HH_Read(Telemetry_t *telemetry);
void IIS2MDCTR_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin);
void IIS2MDCTR_Read(Telemetry_t *telemetry);

uint8_t LSM6DSL_ReadReg(uint8_t reg);
uint8_t LPS22HH_ReadReg(uint8_t reg);
void    LPS22HH_TestRead(float *pressure_hpa);
void ICM45686_TestRead(int16_t *ax, int16_t *ay, int16_t *az);
void IIS2MDCTR_TestRead(int16_t *mx, int16_t *my, int16_t *mz);
uint8_t LSM6DSL_ReadReg(uint8_t reg);
void    LSM6DSL_WriteReg(uint8_t reg, uint8_t val);

