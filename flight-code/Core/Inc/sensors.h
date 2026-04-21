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

// ADXL375BCCZ-RL7 (IMU) Registers
#define ADXL375_DEVID        0x00
#define ADXL375_BW_RATE      0x2C
#define ADXL375_POWER_CTL    0x2D
#define ADXL375_DATA_FORMAT  0x31
#define ADXL375_DATAX0       0x32

// LIS (Mag) Registers
#define LIS3MDLTR_WHO_AM_I    0x0F
#define LIS3MDLTR_CTRL_REG1   0x20
#define LIS3MDLTR_CTRL_REG2   0x21
#define LIS3MDLTR_CTRL_REG3   0x22
#define LIS3MDLTR_CTRL_REG5	  0x24
#define LIS3MDLTR_OUT_X_L     0x28


//BMX055 Registers

//BMP388 Registers
#define BMP388_CHIP_ID     0x00
#define BMP388_STATUS      0x03
#define BMP388_PRESS_DATA  0x04   // 3 bytes
#define BMP388_TEMP_DATA   0x07   // 3 bytes
#define BMP388_PWR_CTRL    0x1B
#define BMP388_OSR         0x1C
#define BMP388_ODR         0x1D
#define BMP388_CONFIG      0x1F
#define BMP388_CALIB_DATA  0x31
#define BMP388_CMD         0x7E

typedef struct {
    float par_t1;
    float par_t2;
    float par_t3;
    float par_p1;
    float par_p2;
    float par_p3;
    float par_p4;
    float par_p5;
    float par_p6;
    float par_p7;
    float par_p8;
    float par_p9;
    float par_p10;
    float par_p11;
    float t_lin;
} BMP388_CalibData;
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

// ADXL375BCCZ-RL7 Functions
void ADXL375_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin);
void ADXL375_Read(Telemetry_t *telemetry);
uint8_t ADXL375_WhoAmI(void);

// LIS3MDLTR Functions
void LIS3MDLTR_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin);
void LIS3MDLTR_Read(Telemetry_t *telemetry);
uint8_t LIS3MDLTR_WhoAmI(void);

//BMX055 Functions

//BMP388 Functions
void BMP388_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin);
void BMP388_Read(Telemetry_t *telemetry);
uint8_t BMP388_WhoAmI(void);

// Test/Debug Functions
uint8_t LPS22HH_WhoAmI(void);
uint8_t LPS22HH_ReadReg(uint8_t reg);
void LPS22HH_WriteReg(uint8_t reg, uint8_t val);

void transform_accel_to_world(Telemetry_t *telemetry);
void deselect_all_spi();

// Bias
void Bias_Init(Bias_t *bias);
void Apply_Bias(Bias_t *bias, Telemetry_t *telemetry);
void Gyro_CalibrateBias(Bias_t* bias, Telemetry_t* telemetry, int num_samples);
void calibrate_accel_bias_stationary(Telemetry_t *telemetry);

void log_mag();
