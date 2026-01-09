#pragma once

#include "telemetry.h"
#include "main.h"

// LSMDSL (IMU 1)
#define LSM6DSL_CTRL1_XL   0x10
#define LSM6DSL_CTRL2_G    0x11
#define LSM6DSL_CTRL3_C    0x12
#define LSM6DSL_OUTX_L_G   0x22

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

void init_sensors();
void read_sensors(Telemetry_t *telemetry);
void read_bmp(Telemetry_t *telemetry);
