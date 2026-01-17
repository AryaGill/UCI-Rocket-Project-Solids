#pragma once

#include "main.h"

// LPS22HHTR (Barometer) Registers
#define LPS22HH_WHO_AM_I    0x0F
#define LPS22HH_CTRL_REG1   0x10
#define LPS22HH_CTRL_REG2   0x11
#define LPS22HH_STATUS_REG  0x27
#define LPS22HH_PRESS_OUT_XL 0x28
#define LPS22HH_TEMP_OUT_L   0x2B

// Function Declarations
void init_sensors(SPI_HandleTypeDef *hspi);
void read_sensors(Telemetry_t *telemetry);

// LPS22HHTR Functions
void LPS22HH_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin);
void LPS22HH_Read(Telemetry_t *telemetry);

// Test/Debug Functions
uint8_t LPS22HH_WhoAmI(void);
uint8_t LPS22HH_ReadReg(uint8_t reg);
void LPS22HH_WriteReg(uint8_t reg, uint8_t val);
