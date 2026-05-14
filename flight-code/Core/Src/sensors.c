#include "sensors.h"
#include "sd_card.h"
#include <stdio.h>
#include <string.h>

GPIO_TypeDef *LPS22HH_port;
uint16_t LPS22HH_pin;
SPI_HandleTypeDef *LPS22HH_hspi;

GPIO_TypeDef *LSM_port;
uint16_t LSM_pin;
SPI_HandleTypeDef *LSM_hspi;

GPIO_TypeDef *LIS_port;
uint16_t LIS_pin;
SPI_HandleTypeDef *LIS_hspi;

GPIO_TypeDef *BMP388_port;
uint16_t BMP388_pin;
SPI_HandleTypeDef *BMP388_hspi;

//BMX
// Accelerometer
GPIO_TypeDef *BMX_ACC_port;
uint16_t BMX_ACC_pin;
SPI_HandleTypeDef *BMX_ACC_hspi;

// Gyroscope
GPIO_TypeDef *BMX_GYRO_port;
uint16_t BMX_GYRO_pin;
SPI_HandleTypeDef *BMX_GYRO_hspi;

// Magnetometer
GPIO_TypeDef *BMX_MAG_port;
uint16_t BMX_MAG_pin;
SPI_HandleTypeDef *BMX_MAG_hspi;

extern Bias_t bias;

volatile uint8_t lps_whoami = 0; // Should be 0xB3 for LPS22HH
volatile uint8_t lsm_whoami = 0; // Should be 0x6A
volatile uint8_t lis_whoami = 0; // Should be 0x3D
volatile uint8_t bmp388_whoami = 0; //Should be 0x50

volatile uint8_t bmx_acc_whoami = 0; // Should be 11111010 or 0xFA
volatile uint8_t bmx_gyro_whoami = 0; // Should be 0x0f
volatile uint8_t bmx_mag_whoami = 0; //Should be 0x32

static BMP388_CalibData bmp388_calib; //calibration struct for bmp388

float max_r;
float max_p;
float max_y;
float min_r;
float min_p;
float min_y;


uint8_t Verify_Sensors(void){
	//Check Barometer
	lps_whoami = LPS22HH_WhoAmI();
	if (lps_whoami != 0xB3){
		return 1;
	}

	//Check LSM6DSL IMU
//	lsm_whoami = LSM6DSL_WhoAmI();
//	if (lsm_whoami != 0x6a){
//		return 1;
//	}
//
//	lis_whoami = LIS3MDLTR_WhoAmI();
//	if (lis_whoami != 0x3D){
//		return 1;
//	}

//	bmp388_whoami = BMP388_WhoAmI();
//	if (bmp388_whoami != 0x50){
//		return 1;
//	}
	//BMX
	bmx_acc_whoami = BMX055_ACC_WhoAmI();
	if (bmx_acc_whoami != 0xFA){
		return 1;
	}
	bmx_gyro_whoami = BMX055_GYRO_WhoAmI();
	if (bmx_gyro_whoami != 0x0f){
		return 1;
	}

//	bmx_mag_whoami = BMX055_MAG_WhoAmI();
//	if (bmx_mag_whoami != 0x32){
//		return 1;
//	}
	return 0;
}

// SPI Helper Functions
static inline void SPI_CS_LOW(GPIO_TypeDef *port, uint16_t pin)
{
    HAL_GPIO_WritePin(port, pin, GPIO_PIN_RESET);
}

static inline void SPI_CS_HIGH(GPIO_TypeDef *port, uint16_t pin)
{
    HAL_GPIO_WritePin(port, pin, GPIO_PIN_SET);
}

static void SPI_Read(SPI_HandleTypeDef *hspi, GPIO_TypeDef *port, uint16_t pin,
                     uint8_t reg, uint8_t *buf, uint8_t len)
{
    reg |= 0x80; // Set read bit
    SPI_CS_LOW(port, pin);
    HAL_SPI_Transmit(hspi, &reg, 1, HAL_MAX_DELAY);
    HAL_SPI_Receive(hspi, buf, len, HAL_MAX_DELAY);
    SPI_CS_HIGH(port, pin);
}

static void SPI_Read_Multi(SPI_HandleTypeDef *hspi, GPIO_TypeDef *port, uint16_t pin,
                     uint8_t reg, uint8_t *buf, uint8_t len)
{
    reg |= 0xC0; // Set read bit and auto increment
    SPI_CS_LOW(port, pin);
    HAL_SPI_Transmit(hspi, &reg, 1, HAL_MAX_DELAY);
    HAL_SPI_Receive(hspi, buf, len, HAL_MAX_DELAY);
    SPI_CS_HIGH(port, pin);
}

static void SPI_Write(SPI_HandleTypeDef *hspi, GPIO_TypeDef *port, uint16_t pin,
                      uint8_t reg, uint8_t val)
{
    uint8_t tx[2] = { reg & 0x7F, val }; // Clear read bit
    SPI_CS_LOW(port, pin);
    HAL_SPI_Transmit(hspi, tx, 2, HAL_MAX_DELAY);
    SPI_CS_HIGH(port, pin);
}

// Sensor Initialization
void init_sensors(SPI_HandleTypeDef *hspi2, SPI_HandleTypeDef *hspi4)
{
    // Force all CS HIGH immediately
    HAL_GPIO_WritePin(Baro_CS_GPIO_Port, Baro_CS_Pin, GPIO_PIN_SET);
    HAL_GPIO_WritePin(IMU_2_CS_GPIO_Port, IMU_2_CS_Pin, GPIO_PIN_SET);
    HAL_GPIO_WritePin(Mag_CS_GPIO_Port, Mag_CS_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(Baro2_CS_GPIO_Port, Baro2_CS_Pin, GPIO_PIN_SET);
    HAL_GPIO_WritePin(BMX_ACCEL_CS_GPIO_Port, BMX_ACCEL_CS_Pin, GPIO_PIN_SET);
    HAL_GPIO_WritePin(BMX_GYRO_CS_GPIO_Port, BMX_GYRO_CS_Pin, GPIO_PIN_SET);
    HAL_GPIO_WritePin(BMX_MAG_CS_GPIO_Port, BMX_MAG_CS_Pin, GPIO_PIN_SET);

    HAL_Delay(100);

    // Initialize Baro 1
    LPS22HH_Init(hspi2, Baro_CS_GPIO_Port, Baro_CS_Pin);
    HAL_Delay(20);

    // Initialize BMX 1
    BMX055_Init(hspi2, BMX_ACCEL_CS_GPIO_Port, BMX_ACCEL_CS_Pin, BMX_GYRO_CS_GPIO_Port, BMX_GYRO_CS_Pin, BMX_MAG_CS_GPIO_Port, BMX_MAG_CS_Pin);

    Bias_Init(&bias);

//    // Initialize LIS
//    LIS3MDLTR_Init(hspi4, Mag_CS_GPIO_Port, Mag_CS_Pin);
//    HAL_Delay(20);
}

// Sensor Reading
void read_sensors(Telemetry_t *telemetry)
{
    LPS22HH_Read(telemetry);
//    LIS3MDLTR_Read(telemetry);
    BMX055_Read(telemetry);

//    Apply_Bias(&bias, telemetry);

//    transform_accel_to_world(telemetry);

    telemetry->time = HAL_GetTick();

    // Comment out for flight
//    calibrate_accel_bias_stationary(telemetry);
}

// Calculate altitude from pressure (standard atmosphere model)

float Calculate_Altitude(float pressure_hPa)
{
	const float sea_level_pressure = 1013.25f; // hPa at sea level

	// Barometric formula: h = 44330 * (1 - (P/P0)^(1/5.255))
	float ratio = pressure_hPa / sea_level_pressure;
	float altitude = 44330.0f * (1.0f - powf(ratio, 0.1903f));

	return altitude;
}

// LPS22HHTR Functions
void LPS22HH_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin)
{
    LPS22HH_port = cs_port;
    LPS22HH_pin = cs_pin;
    LPS22HH_hspi = hspi;

    HAL_Delay(20); // Wait for sensor power-up

    lps_whoami = LPS22HH_WhoAmI();

    // Software reset
    SPI_Write(hspi, cs_port, cs_pin, LPS22HH_CTRL_REG2, 0x04);
    HAL_Delay(10);

    // Configure: 75 Hz ODR, continuous mode, BDU enabled
    // CTRL_REG1: ODR=75Hz (0101), AVG=512 (11), EN_LPFP=0, BDU=1
    SPI_Write(hspi, cs_port, cs_pin, LPS22HH_CTRL_REG1, 0x5C);
    HAL_Delay(10);
}

void LPS22HH_Read(Telemetry_t *telemetry)
{
	uint8_t buf[5];
	int32_t raw_p;
    int16_t raw_t;

    // Read pressure (3 bytes) and temperature (2 bytes) - 5 bytes total
    SPI_Read(LPS22HH_hspi, LPS22HH_port, LPS22HH_pin, LPS22HH_PRESS_OUT_XL, buf, 5);

    // Pressure is 24-bit, little-endian
    raw_p = (int32_t)(((uint32_t)buf[2] << 16) | ((uint32_t)buf[1] << 8) | buf[0]);

    // Sign extend from 24-bit to 32-bit
    if (raw_p & 0x800000) {
        raw_p |= 0xFF000000;
    }

    // Temperature is 16-bit, little-endian
    raw_t = (int16_t)(((uint16_t)buf[4] << 8) | buf[3]);

    // Convert to metric units
    telemetry->pressure = raw_p / 4096.0f;      // hPa (mbar)
    telemetry->temperature = raw_t / 100.0f;     // °C
    telemetry->altitude = Calculate_Altitude(telemetry->pressure);
}

uint8_t LPS22HH_WhoAmI(void)
{
    uint8_t id = 0;
    SPI_Read(LPS22HH_hspi,
             LPS22HH_port,
             LPS22HH_pin,
             LPS22HH_WHO_AM_I,
             &id,
             1);
    return id;
}

uint8_t LPS22HH_ReadReg(uint8_t reg)
{
    uint8_t val = 0;
    SPI_Read(LPS22HH_hspi,
             LPS22HH_port,
             LPS22HH_pin,
             reg,
             &val,
             1);
    return val;
}

void LPS22HH_WriteReg(uint8_t reg, uint8_t val)
{
    SPI_Write(LPS22HH_hspi,
              LPS22HH_port,
              LPS22HH_pin,
              reg,
              val);
}

void LSM6DSL_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin) {
    LSM_port = cs_port;
    LSM_pin = cs_pin;
    LSM_hspi = hspi;

    // 1. Software Reset
    SPI_Write(hspi, cs_port, cs_pin, LSM6DSL_CTRL3_C, 0x01);
    HAL_Delay(50);

    // 2. Enable Block Data Update (BDU) and Auto-Increment
    // CTRL3_C: BDU=1 (0x40), IF_INC=1 (0x04) -> 0x44
    SPI_Write(hspi, cs_port, cs_pin, LSM6DSL_CTRL3_C, 0x44);

    // 3. Configure Accelerometer: 104Hz, +/- 4g
    // CTRL1_XL: 0100 (104Hz), 01 (16g) -> 0x48
    SPI_Write(hspi, cs_port, cs_pin, LSM6DSL_CTRL1_XL, 0x44);

    // 4. Configure Gyroscope: 104Hz, 2000 dps
    // CTRL2_G: 0100 (104Hz), 11 (2000dps) -> 0x4C
    SPI_Write(hspi, cs_port, cs_pin, LSM6DSL_CTRL2_G, 0x4C);
}

void LSM6DSL_Read(Telemetry_t *telemetry) {
    uint8_t buf[12];
    // Read 12 bytes starting from Gyro X_L all the way through Accel Z_H
    SPI_Read(LSM_hspi, LSM_port, LSM_pin, LSM6DSL_OUTX_L_G, buf, 12);

    // Raw values (Little Endian)
    int16_t gx = (int16_t)((buf[1] << 8) | buf[0]);
    int16_t gy = (int16_t)((buf[3] << 8) | buf[2]);
    int16_t gz = (int16_t)((buf[5] << 8) | buf[4]);
    int16_t ax = (int16_t)((buf[7] << 8) | buf[6]);
    int16_t ay = (int16_t)((buf[9] << 8) | buf[8]);
    int16_t az = (int16_t)((buf[11] << 8) | buf[10]);

    // Accelerometer (±16 g)
	// 0.488 mg/LSB → 0.000488 g/LSB
	// Convert to m/s²: * 9.80665
	const float ACCEL_SCALE = 0.000488f * 9.80665f;  // ≈ 0.00479

	telemetry->lsm_accel_r = ay * ACCEL_SCALE;
	telemetry->lsm_accel_p = ax * ACCEL_SCALE;
	telemetry->lsm_accel_y = -az * ACCEL_SCALE;

	// Gyroscope (2000 dps)
	// 70 mdps/LSB = 0.07 dps/LSB
	// Convert to rad/s: * (π / 180)
	const float GYRO_SCALE = 0.07f * (3.14159265359f / 180.0f); // ≈ 0.00122173

	telemetry->lsm_gyro_r = gy * GYRO_SCALE;
	telemetry->lsm_gyro_p = gx * GYRO_SCALE;
	telemetry->lsm_gyro_y = -gz * GYRO_SCALE;
}

uint8_t LSM6DSL_WhoAmI(void) {
    uint8_t id = 0;
    SPI_Read(LSM_hspi, LSM_port, LSM_pin, LSM6DSL_WHO_AM_I, &id, 1);
    return id;
}

void LIS3MDLTR_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin) {
    LIS_port = cs_port;
    LIS_pin = cs_pin;
    LIS_hspi = hspi;

    // 1. Reset the device (CTRL_REG2: soft reset)
    SPI_Write(hspi, cs_port, cs_pin, LIS3MDLTR_CTRL_REG2, 0x0C); // soft reset + reboot
    HAL_Delay(50);

    lis_whoami = LIS3MDLTR_WhoAmI();

    // 2. Enable Block Data Update (BDU) and continuous mode
    // CTRL_REG1: Temp sensor off, Ultra-high-performance XY, 80Hz ODR
    SPI_Write(hspi, cs_port, cs_pin, LIS3MDLTR_CTRL_REG1, 0x70); // 0111 0000

    // 3. Set full-scale to ±4 gauss, ultra-high-performance Z
    SPI_Write(hspi, cs_port, cs_pin, LIS3MDLTR_CTRL_REG2, 0x00); // FS = ±4 gauss

    // 4. Enable continuous-conversion mode
    SPI_Write(hspi, cs_port, cs_pin, LIS3MDLTR_CTRL_REG3, 0x00); // Continuous-conversion
}

void LIS3MDLTR_Read(Telemetry_t *telemetry) {
	uint8_t buf[6];
    // Read 6 bytes starting from OUT_X_L
	SPI_Read_Multi(LIS_hspi, LIS_port, LIS_pin, LIS3MDLTR_OUT_X_L, buf, 6);

    // Raw values (Little Endian)
    int16_t mx = (int16_t)((buf[1] << 8) | buf[0]);
    int16_t my = (int16_t)((buf[3] << 8) | buf[2]);
    int16_t mz = (int16_t)((buf[5] << 8) | buf[4]);

    // Conversion to microteslas
    // ±4 gauss full-scale = 0.14 mG/LSB = 0.014 µT/LSB * 100? Actually LIS3MDLTR FS=4G => 0.14 mG/LSB
    // Let's compute factor: 1 G = 100 µT, so 0.14 mG = 0.014 µT
    float factor = 0.014f;

    telemetry->mag_r = -mx * factor;
    telemetry->mag_p = my * factor;
    telemetry->mag_y = -mz * factor;
}

uint8_t LIS3MDLTR_WhoAmI(void) {
    uint8_t id = 0;
    SPI_Read(LIS_hspi, LIS_port, LIS_pin, LIS3MDLTR_WHO_AM_I, &id, 1);
    return id;
}


void BMX055_Init(SPI_HandleTypeDef *hspi,
                 GPIO_TypeDef *acc_port, uint16_t acc_pin,
                 GPIO_TypeDef *gyro_port, uint16_t gyro_pin,
                 GPIO_TypeDef *mag_port, uint16_t mag_pin)
{
    BMX_ACC_hspi = hspi;
    BMX_ACC_port = acc_port;
    BMX_ACC_pin = acc_pin;

    BMX_GYRO_hspi = hspi;
    BMX_GYRO_port = gyro_port;
    BMX_GYRO_pin = gyro_pin;

    BMX_MAG_hspi = hspi;
    BMX_MAG_port = mag_port;
    BMX_MAG_pin = mag_pin;

    HAL_Delay(10);

    // ---------- ACCEL ----------
    SPI_Write(BMX_ACC_hspi, BMX_ACC_port, BMX_ACC_pin,
              BMX055_ACC_SOFTRESET, BMX055_ACC_SOFTRESET_CMD);
    HAL_Delay(10);

    bmx_acc_whoami = BMX055_ACC_WhoAmI();


    SPI_Write(BMX_ACC_hspi, BMX_ACC_port, BMX_ACC_pin,
              BMX055_ACC_PMU_RANGE, BMX055_ACC_RANGE_16G);

    SPI_Write(BMX_ACC_hspi, BMX_ACC_port, BMX_ACC_pin,
              BMX055_ACC_PMU_BW, BMX055_ACC_BW_125HZ);

    // ---------- GYRO ----------
    SPI_Write(BMX_GYRO_hspi, BMX_GYRO_port, BMX_GYRO_pin,
              BMX055_GYRO_SOFTRESET, BMX055_GYRO_SOFTRESET_CMD);
    HAL_Delay(10);

    bmx_gyro_whoami = BMX055_GYRO_WhoAmI();

    SPI_Write(BMX_GYRO_hspi, BMX_GYRO_port, BMX_GYRO_pin,
              BMX055_GYRO_RANGE, BMX055_GYRO_RANGE_500DPS);

    SPI_Write(BMX_GYRO_hspi, BMX_GYRO_port, BMX_GYRO_pin,
              BMX055_GYRO_BW, BMX055_GYRO_BW_200HZ);

//    // ---------- MAG ----------
    SPI_Write(BMX_MAG_hspi, BMX_MAG_port, BMX_MAG_pin,
              BMX055_MAG_POWER_CTRL, BMX055_MAG_POWER_ON);
    HAL_Delay(10);

    bmx_mag_whoami = BMX055_MAG_WhoAmI();

    SPI_Write(BMX_MAG_hspi, BMX_MAG_port, BMX_MAG_pin,
              BMX055_MAG_OP_MODE, BMX055_MAG_NORMAL_MODE);

    SPI_Write(BMX_MAG_hspi, BMX_MAG_port, BMX_MAG_pin,
              BMX055_MAG_REP_XY, BMX055_MAG_REPXY_DEFAULT);

    SPI_Write(BMX_MAG_hspi, BMX_MAG_port, BMX_MAG_pin,
              BMX055_MAG_REP_Z, BMX055_MAG_REPZ_DEFAULT);

    HAL_Delay(10);
}

uint8_t BMX055_ACC_WhoAmI(void)
{
    uint8_t id;
    SPI_Read(BMX_ACC_hspi, BMX_ACC_port, BMX_ACC_pin,
             BMX055_ACC_CHIP_ID, &id, 1);
    return id;
}

uint8_t BMX055_GYRO_WhoAmI(void)
{
    uint8_t id;
    SPI_Read(BMX_GYRO_hspi, BMX_GYRO_port, BMX_GYRO_pin,
             BMX055_GYRO_CHIP_ID, &id, 1);
    return id;
}

uint8_t BMX055_MAG_WhoAmI(void)
{
    uint8_t id;
    SPI_Read(BMX_MAG_hspi, BMX_MAG_port, BMX_MAG_pin,
             BMX055_MAG_CHIP_ID, &id, 1);
    return id;
}

void BMX055_Read_Accel(Telemetry_t *t)
{
    uint8_t buf[6];
    SPI_Read(BMX_ACC_hspi, BMX_ACC_port, BMX_ACC_pin,
             BMX055_ACC_X_LSB, buf, 6);

    int16_t ax = ((int16_t)((int16_t)buf[1] << 8 | buf[0])) >> 4;
    int16_t ay = ((int16_t)((int16_t)buf[3] << 8 | buf[2])) >> 4;
    int16_t az = ((int16_t)((int16_t)buf[5] << 8 | buf[4])) >> 4;

    const float scale = 9.81f / 128.0f;

    t->bmx_accel_r = -ay * scale;
    t->bmx_accel_p = ax * scale;
    t->bmx_accel_y = az * scale;
}

void BMX055_Read_Gyro(Telemetry_t *t)
{
    uint8_t buf[6];
    SPI_Read(BMX_GYRO_hspi, BMX_GYRO_port, BMX_GYRO_pin,
             BMX055_GYRO_RATE_X_LSB, buf, 6);

    int16_t gx = (int16_t)((int16_t)buf[1] << 8 | buf[0]);
    int16_t gy = (int16_t)((int16_t)buf[3] << 8 | buf[2]);
    int16_t gz = (int16_t)((int16_t)buf[5] << 8 | buf[4]);

    const float scale = (1.0f / 65.5f) * (M_PI / 180.0f); //radians

    t->bmx_gyro_r = -gy * scale;
    t->bmx_gyro_p = gx * scale;
    t->bmx_gyro_y = gz * scale;
}

void BMX055_Read_Mag(Telemetry_t *t)
{
    uint8_t buf[8];

    SPI_Read(BMX_MAG_hspi, BMX_MAG_port, BMX_MAG_pin,
             BMX055_MAG_DATA_X_LSB, buf, 8);

    // Cast buf[n] to int16_t BEFORE shifting to preserve sign extension
    int16_t mx = ((int16_t)((int16_t)buf[1] << 8 | buf[0])) >> 3;
    int16_t my = ((int16_t)((int16_t)buf[3] << 8 | buf[2])) >> 3;
    int16_t mz = ((int16_t)((int16_t)buf[5] << 8 | buf[4])) >> 1;

    t->bmx_mag_r = -mx;
    t->bmx_mag_p = my;
    t->bmx_mag_y = mz;
}

void BMX055_Read(Telemetry_t *t)
{
    BMX055_Read_Accel(t);
    BMX055_Read_Gyro(t);
    BMX055_Read_Mag(t);
}

void calibrate_mag(Telemetry_t *telemetry){
//	LIS3MDLTR_Read(telemetry);
//
//	// Initialize max and min
//	if (max_r == 0) max_r = telemetry->mag_r;
//	if (max_p == 0) max_p = telemetry->mag_p;
//	if (max_y == 0) max_y = telemetry->mag_y;
//	if (min_r == 0) min_r = telemetry->mag_r;
//	if (min_y == 0) min_y = telemetry->mag_y;
//	if (min_p == 0) min_p = telemetry->mag_p;
//
//	// Update max and min
//	if (telemetry->mag_r > max_r) max_r = telemetry->mag_r;
//	if (telemetry->mag_p > max_p) max_p = telemetry->mag_p;
//	if (telemetry->mag_y > max_y) max_y = telemetry->mag_y;
//	if (telemetry->mag_r < min_r) min_r = telemetry->mag_r;
//	if (telemetry->mag_p < min_p) min_p = telemetry->mag_p;
//	if (telemetry->mag_y < min_y) min_y = telemetry->mag_y;
//
//	float offset_r = (max_r + min_r) / 2.0f;
//	float offset_p = (max_p + min_p) / 2.0f;
//	float offset_y = (max_y + min_y) / 2.0f;
//
//	float radius_r = (max_r - min_r) / 2.0f;
//	float radius_p = (max_p - min_p) / 2.0f;
//	float radius_y = (max_y - min_y) / 2.0f;
//
//	float avg_radius = (radius_r + radius_p + radius_y) / 3.0f;
//
//	float scale_r = avg_radius / radius_r;
//	float scale_p = avg_radius / radius_p;
//	float scale_y = avg_radius / radius_y;
//
//	write_mag("MAG_CALIB.csv", offset_r, offset_p, offset_y, scale_r, scale_p, scale_y);
}

void transform_accel_to_world(Telemetry_t *telemetry) {
  // Average IMUs (body frame)
  float ax = telemetry->bmx_accel_p;
  float ay = telemetry->bmx_accel_y;
  float az = telemetry->bmx_accel_r;

  // Quaternion (w, x, y, z)
  float qw = telemetry->q0;
  float qx = telemetry->q1;
  float qy = telemetry->q2;
  float qz = telemetry->q3;

  // Rotation matrix (body → world)
  float R11 = 1.0f - 2.0f*(qy*qy + qz*qz);
  float R12 = 2.0f*(qx*qy - qz*qw);
  float R13 = 2.0f*(qx*qz + qy*qw);

  float R21 = 2.0f*(qx*qy + qz*qw);
  float R22 = 1.0f - 2.0f*(qx*qx + qz*qz);
  float R23 = 2.0f*(qy*qz - qx*qw);

  float R31 = 2.0f*(qx*qz - qy*qw);
  float R32 = 2.0f*(qy*qz + qx*qw);
  float R33 = 1.0f - 2.0f*(qx*qx + qy*qy);

  // Rotate acceleration into world frame
  telemetry->accel_world_x = R11*ax + R12*ay + R13*az;
  telemetry->accel_world_y = R21*ax + R22*ay + R23*az;
  telemetry->accel_world_z = R31*ax + R32*ay + R33*az - 9.81f;
}

void deselect_all_spi(){
	HAL_GPIO_WritePin(RF_CS_GPIO_Port, RF_CS_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(IMU_2_CS_GPIO_Port, IMU_2_CS_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(Mag_CS_GPIO_Port, Mag_CS_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(Baro_CS_GPIO_Port, Baro_CS_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(Flash_CS_GPIO_Port, Flash_CS_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(SD_CS_GPIO_Port, SD_CS_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(Baro2_CS_GPIO_Port, Baro2_CS_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(BMX_ACCEL_CS_GPIO_Port, BMX_ACCEL_CS_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(BMX_GYRO_CS_GPIO_Port, BMX_GYRO_CS_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(BMX_MAG_CS_GPIO_Port, BMX_MAG_CS_Pin, GPIO_PIN_SET);
}

void Gyro_CalibrateBias(Bias_t* bias, Telemetry_t* telemetry, int num_samples){

	// Get avg gyro values give num_samples.
	// WARNING! MUST BE DONE WHILE STATIONARY


	// Init calibration values at 0
	float sum_gr = 0.0f;
	float sum_gp = 0.0f;
	float sum_gy = 0.0f;

	// Sum all values for num_smaples
	for (int i = 0; i < num_samples; i++){
		read_sensors(telemetry);
		// converted from lsm -> bmx
		sum_gr += telemetry->bmx_gyro_r + bias->bmx_gyro_r_bias;
		sum_gp += telemetry->bmx_gyro_p + bias->bmx_gyro_p_bias;
		sum_gy += telemetry->bmx_gyro_y + bias->bmx_gyro_y_bias;

		HAL_Delay(2);
	}

	// Divide by total samples to get average
	bias->bmx_gyro_r_bias = sum_gr / num_samples;
	bias->bmx_gyro_p_bias = sum_gp / num_samples;
	bias->bmx_gyro_y_bias = sum_gy / num_samples;
}

void Bias_Init(Bias_t *bias)
{
	bias->bmx_gyro_r_bias = 0.0;
	bias->bmx_gyro_p_bias = 0.0;
	bias->bmx_gyro_y_bias = 0.0f;

	bias->bmx_accel_r_bias = 0.0f;
	bias->bmx_accel_p_bias = 0.0f;
	bias->bmx_accel_y_bias = 0.0f;
	bias->bias_count = 0;

    bias->mag_r_bias = 0.0f; // 25.83, 25.305, 26.313
    bias->mag_p_bias = 0.0f; // 3.934, 1.547, 1.106
    bias->mag_y_bias = 0.0f; // 22.043, 21.798, 18.046

    bias->mag_r_scale = 1.0;
    bias->mag_p_scale = 1.0;
    bias->mag_y_scale = 1.0;
}

void Bias_Calculate(Bias_t *bias, Telemetry_t *t, int num_samples){

	Gyro_CalibrateBias(bias, t, num_samples);


	bias->bias_count = num_samples;

//    bias->bias_count += 1.0f;
//
//    float n = bias->bias_count;
//

//
//    bias->lsm_accel_r_bias += ((t->lsm_accel_r - 9.81) - bias->lsm_accel_r_bias) / n;
//    bias->lsm_accel_p_bias += (t->lsm_accel_p - bias->lsm_accel_p_bias) / n;
//    bias->lsm_accel_y_bias += (t->lsm_accel_y - bias->lsm_accel_y_bias) / n;
}

void Apply_Bias(Bias_t *bias, Telemetry_t *t)
{

	t->bmx_gyro_r -= bias->bmx_gyro_r_bias;
	t->bmx_gyro_p -= bias->bmx_gyro_p_bias;
	t->bmx_gyro_y -= bias->bmx_gyro_y_bias;

	t->bmx_accel_r -= bias->bmx_accel_r_bias;
	t->bmx_accel_p -= bias->bmx_accel_p_bias;
	t->bmx_accel_y -= bias->bmx_accel_y_bias;


    t->bmx_mag_r -= bias->mag_r_bias;
    t->bmx_mag_p -= bias->mag_p_bias;
    t->bmx_mag_y -= bias->mag_y_bias;

    t->bmx_mag_r *= bias->mag_r_scale;
    t->bmx_mag_p *= bias->mag_p_scale;
    t->bmx_mag_y *= bias->mag_y_scale;
}

// How to calibrate accel: run in debugger. Hold still.
// Place the board in 6 different orientations, letting gravity (1g) act on each axis positively and negatively.
// +X, -X, +Y, -Y, +Z, -Z facing down.
// For each position, record the average output of all three axes.
// For each axis, the offset is: Offset = (Value_+1g + Value_-1g) / 2
// The scale factor is: Scale = (Value_+1g - Value_-1g) / 2 (Theoretically, this should be 1g, but you can use it to correct minor gain errors).
float bmx_accel_r_bias_instance;
float bmx_accel_p_bias_instance;
float bmx_accel_y_bias_instance;
uint32_t prev_time_accel_bias;

void calibrate_accel_bias_stationary(Telemetry_t *telemetry)
{
	// Average IMUs (body frame)
	float ax = telemetry->bmx_accel_p;
	float ay = telemetry->bmx_accel_y;
	float az = telemetry->bmx_accel_r;

	// Quaternion (w, x, y, z)
	float qw = telemetry->q0;
	float qx = telemetry->q1;
	float qy = telemetry->q2;
	float qz = telemetry->q3;

    // Rotation matrix (body -> world)
    float R11 = 1.0f - 2.0f*(qy*qy + qz*qz);
    float R12 = 2.0f*(qx*qy - qz*qw);
    float R13 = 2.0f*(qx*qz + qy*qw);

    float R21 = 2.0f*(qx*qy + qz*qw);
    float R22 = 1.0f - 2.0f*(qx*qx + qz*qz);
    float R23 = 2.0f*(qy*qz - qx*qw);

    float R31 = 2.0f*(qx*qz - qy*qw);
    float R32 = 2.0f*(qy*qz + qx*qw);
    float R33 = 1.0f - 2.0f*(qx*qx + qy*qy);

    // Step 3: measured accel in world frame
    float ax_w = R11*ax + R12*ay + R13*az;
    float ay_w = R21*ax + R22*ay + R23*az;
    float az_w = R31*ax + R32*ay + R33*az;

    // Step 4: subtract gravity (expected world accel = [0,0,9.81])
    float err_wx = ax_w;
    float err_wy = ay_w;
    float err_wz = az_w - 9.81f;

    // Step 5: rotate error back into body frame (R^T)
    float bias_x = R11*err_wx + R21*err_wy + R31*err_wz;
    float bias_y = R12*err_wx + R22*err_wy + R32*err_wz;
    float bias_z = R13*err_wx + R23*err_wy + R33*err_wz;

    float TAU = 1;
    uint32_t now = micros();
    float dt = (now - prev_time_accel_bias) * 1e-6;
    prev_time_accel_bias = now;
    float alpha = TAU / (TAU + dt);
    bmx_accel_p_bias_instance = bmx_accel_p_bias_instance * alpha + bias_x * (alpha - 1);
    bmx_accel_y_bias_instance = bmx_accel_y_bias_instance * alpha + bias_y * (alpha - 1);
    bmx_accel_r_bias_instance = bmx_accel_r_bias_instance * alpha + bias_z * (alpha - 1);
}
