#pragma once

#include "main.h"

// LPS22HHTR (Barometer) Registers
#define LPS22HH_WHO_AM_I    0x0F
#define LPS22HH_CTRL_REG1   0x10
#define LPS22HH_CTRL_REG2   0x11
#define LPS22HH_STATUS_REG  0x27
#define LPS22HH_PRESS_OUT_XL 0x28
#define LPS22HH_TEMP_OUT_L   0x2B

// LSM6DSL (IMU) Registers
#define LSM6DSL_WHO_AM_I            0x0F
#define LSM6DSL_CTRL1_XL            0x10
#define LSM6DSL_CTRL2_G             0x11
#define LSM6DSL_CTRL3_C             0x12
#define LSM6DSL_STATUS_REG          0x1E
#define LSM6DSL_OUTX_L_G            0x22
#define LSM6DSL_OUTX_L_XL           0x28

// ICM45 (IMU2) Registers
#define ICM45_REG_WHO_AM_I          0x72  // Should return 0xE9
#define ICM45_REG_REG_BANK_SEL      0x76
#define ICM45_REG_INTF_CONFIG1      0x4D
#define ICM45_REG_PWR_MGMT0         0x10
#define ICM45_REG_ACCEL_DATA_X1     0x03
#define ICM45_REG_GYRO_DATA_X1      0x09

// Function Declarations
uint8_t Verify_Sensors(void);

void init_sensors(SPI_HandleTypeDef *hspi);
void read_sensors(Telemetry_t *telemetry);

// LPS22HHTR Functions
float Calculate_Altitude(float pressure_hPa);
void LPS22HH_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin);
void LPS22HH_Read(Telemetry_t *telemetry);
void set_start_alt(Telemetry_t *telemetry);

// LSM6DSL Functions
void LSM6DSL_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin);
void LSM6DSL_Read(Telemetry_t *telemetry);
uint8_t LSM6DSL_WhoAmI(void);

// ICM45 Functions
void ICM45686_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin);
void ICM45686_Read(Telemetry_t *telemetry);
uint8_t ICM45686_WhoAmI(void);

// IIS2MDCTR Functions
void IIS2MDCTR_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin, GPIO_TypeDef *sdio_port, uint16_t sdio_pin);
void IIS2MDCTR_Read(Telemetry_t *telemetry);

// Test/Debug Functions
uint8_t LPS22HH_WhoAmI(void);
uint8_t LPS22HH_ReadReg(uint8_t reg);
void LPS22HH_WriteReg(uint8_t reg, uint8_t val);
