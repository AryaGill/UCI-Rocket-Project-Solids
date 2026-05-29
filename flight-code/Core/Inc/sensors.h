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

//BMX055 Registers
//accel
#define BMX055_ACC_CHIP_ID        0x00
#define BMX055_ACC_X_LSB          0x02
#define BMX055_ACC_PMU_RANGE      0x0F
#define BMX055_ACC_PMU_BW         0x10
#define BMX055_ACC_SOFTRESET      0x14

#define BMX055_ACC_RANGE_16G      0x0c
#define BMX055_ACC_BW_125HZ       0x0C
#define BMX055_ACC_SOFTRESET_CMD  0xB6
//gyro
#define BMX055_GYRO_CHIP_ID       0x00
#define BMX055_GYRO_RATE_X_LSB    0x02
#define BMX055_GYRO_RANGE         0x0F
#define BMX055_GYRO_BW            0x10
#define BMX055_GYRO_SOFTRESET     0x14
#define BMX055_GYRO_RANGE_500DPS  0x02
#define BMX055_GYRO_BW_200HZ      0x06
#define BMX055_GYRO_SOFTRESET_CMD 0xB6

//mag
#define BMX055_MAG_CHIP_ID        0x40
#define BMX055_MAG_DATA_X_LSB     0x42
#define BMX055_MAG_DATA_READY     0x48
#define BMX055_MAG_POWER_CTRL     0x4B
#define BMX055_MAG_OP_MODE        0x4C
#define BMX055_MAG_REP_XY         0x51
#define BMX055_MAG_REP_Z          0x52
#define BMX055_MAG_POWER_ON       0x01
#define BMX055_MAG_NORMAL_MODE    0x00
#define BMX055_MAG_REPXY_DEFAULT  0x04
#define BMX055_MAG_REPZ_DEFAULT   0x0F

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

void init_sensors(SPI_HandleTypeDef *hspi2, SPI_HandleTypeDef *hspi4);
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
void BMX055_Init(SPI_HandleTypeDef *hspi,
                 GPIO_TypeDef *acc_port, uint16_t acc_pin,
                 GPIO_TypeDef *gyro_port, uint16_t gyro_pin,
                 GPIO_TypeDef *mag_port, uint16_t mag_pin);

void BMX055_Read(Telemetry_t *telemetry);

void BMX055_Read_Accel(Telemetry_t *telemetry);
void BMX055_Read_Gyro(Telemetry_t *telemetry);
void BMX055_Read_Mag(Telemetry_t *telemetry);

uint8_t BMX055_ACC_WhoAmI(void);
uint8_t BMX055_GYRO_WhoAmI(void);
uint8_t BMX055_MAG_WhoAmI(void);

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

// Bias calculations and application
void Bias_Init(Bias_t *bias);
void Apply_Bias(Bias_t *bias, Telemetry_t *telemetry);
void Gyro_CalibrateBias(Bias_t* bias, Telemetry_t* telemetry, int num_samples);
void calibrate_accel_bias_stationary(Telemetry_t *telemetry);

void log_mag();
