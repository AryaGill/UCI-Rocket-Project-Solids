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

GPIO_TypeDef *ADXL_port;
uint16_t ADXL_pin;
SPI_HandleTypeDef *ADXL_hspi;

GPIO_TypeDef *LIS_port;
uint16_t LIS_pin;
SPI_HandleTypeDef *LIS_hspi;

GPIO_TypeDef *BMP388_port;
uint16_t BMP388_pin;
SPI_HandleTypeDef *BMP388_hspi;


extern Bias_t bias;

volatile uint8_t lps_whoami = 0; // Should be 0xB3 for LPS22HH
volatile uint8_t lsm_whoami = 0; // Should be 0x6A
volatile uint8_t adxl_whoami = 0; // Should be 0xE5
volatile uint8_t lis_whoami = 0; // Should be 0x3D
volatile uint8_t bmp388_whoami = 0; //Should be 0x50

static BMP388_CalibData bmp388_calib; //calibration struct for bmp388

float max_r;
float max_p;
float max_y;
float min_r;
float min_p;
float min_y;


uint8_t Verify_Sensors(void){
	// Check Barometer
	lps_whoami = LPS22HH_WhoAmI();
	if (lps_whoami != 0xB3){
		return 1;
	}

	// Check LSM6DSL IMU
	lsm_whoami = LSM6DSL_WhoAmI();
	if (lsm_whoami != 0x6a){
		return 1;
	}

	adxl_whoami = ADXL375_WhoAmI();
	if (adxl_whoami != 0xE5){
		return 1;
	}

	lis_whoami = LIS3MDLTR_WhoAmI();
	if (lis_whoami != 0x3D){
		return 1;
	}

	bmp388_whoami = BMP388_WhoAmI();
		if (bmp388_whoami != 0x50){
			return 1;
		}

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
void init_sensors(SPI_HandleTypeDef *hspi)
{
    // Force all CS HIGH immediately
    HAL_GPIO_WritePin(Baro_CS_GPIO_Port, Baro_CS_Pin, GPIO_PIN_SET);
    HAL_GPIO_WritePin(Baro2_CS_GPIO_Port, Baro2_CS_Pin, GPIO_PIN_SET);
    HAL_GPIO_WritePin(IMU_CS_GPIO_Port, IMU_CS_Pin, GPIO_PIN_SET);
    HAL_GPIO_WritePin(IMU_2_CS_GPIO_Port, IMU_2_CS_Pin, GPIO_PIN_SET);
    HAL_GPIO_WritePin(Mag_CS_GPIO_Port, Mag_CS_Pin, GPIO_PIN_SET);

    HAL_Delay(100);

    // Initialize Baro
    LPS22HH_Init(hspi, Baro_CS_GPIO_Port, Baro_CS_Pin);
    HAL_Delay(20);

    //Initialize BMP baro
    BMP388_Init(hspi, Baro2_CS_GPIO_Port, Baro2_CS_Pin);
    HAL_Delay(20);

    // Initialize LSM
    LSM6DSL_Init(hspi, IMU_2_CS_GPIO_Port, IMU_2_CS_Pin);
    HAL_Delay(20);

    // Initialize ADXL
    ADXL375_Init(hspi, IMU_CS_GPIO_Port, IMU_CS_Pin);
    HAL_Delay(20);

    // Initialize LIS
    LIS3MDLTR_Init(hspi, Mag_CS_GPIO_Port, Mag_CS_Pin);
    HAL_Delay(20);


    Bias_Init(&bias);
}

// Sensor Reading
void read_sensors(Telemetry_t *telemetry)
{
    LPS22HH_Read(telemetry);
    LSM6DSL_Read(telemetry);
    ADXL375_Read(telemetry);
    LIS3MDLTR_Read(telemetry);
    BMP388_Read(telemetry);

    Apply_Bias(&bias, telemetry);

    transform_accel_to_world(telemetry);

    telemetry->time = HAL_GetTick();

    // Comment out for flight
    calibrate_accel_bias_stationary(telemetry);
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
    raw_p = (int32_t)((buf[2] << 16) | (buf[1] << 8) | buf[0]);

    // Sign extend from 24-bit to 32-bit
    if (raw_p & 0x800000) {
        raw_p |= 0xFF000000;
    }

    // Temperature is 16-bit, little-endian
    raw_t = (int16_t)((buf[4] << 8) | buf[3]);

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

/**
 * Switch SPI1 to Mode 3 (CPOL=1, CPHA=1)
 * Saves previous CPOL/CPHA settings
 */
static inline void SPI_SwitchToMode3(void)
{
//	// Disable SPI before changing mode
//	CLEAR_BIT(SPI1->CR1, SPI_CR1_SPE);
//
//	// Save current CPOL/CPHA (in CFG2 register)
//	spi_saved_mode = SPI1->CFG2 & (SPI_CFG2_CPOL | SPI_CFG2_CPHA);
//
//    // Set Mode 3 (CPOL=1, CPHA=1)
//    SET_BIT(SPI1->CFG2, SPI_CFG2_CPOL | SPI_CFG2_CPHA);
//
//    // Re-enable SPI
//    SET_BIT(SPI1->CR1, SPI_CR1_SPE);

	HAL_SPI_DeInit(ADXL_hspi); // Disable SPI and clean up
	ADXL_hspi->Init.CLKPolarity = SPI_POLARITY_HIGH; // CPOL 1
	ADXL_hspi->Init.CLKPhase = SPI_PHASE_2EDGE;      // CPHA 1
	HAL_SPI_Init(ADXL_hspi);   // Re-initialize with new settings
}

/**
 * Restore previous SPI1 CPOL/CPHA settings
 */
static inline void SPI_RestoreMode(void)
{
//    // Disable SPI before restoring
//    CLEAR_BIT(SPI1->CR1, SPI_CR1_SPE);
//
//    // Restore saved CPOL/CPHA bits
//    MODIFY_REG(SPI1->CFG2,
//               SPI_CFG2_CPOL | SPI_CFG2_CPHA,
//               spi_saved_mode);
//
//    // Re-enable SPI
//    SET_BIT(SPI1->CR1, SPI_CR1_SPE);

	HAL_SPI_DeInit(ADXL_hspi); // Disable SPI and clean up
	ADXL_hspi->Init.CLKPolarity = SPI_POLARITY_LOW;
	ADXL_hspi->Init.CLKPhase = SPI_PHASE_1EDGE;
	HAL_SPI_Init(ADXL_hspi);   // Re-initialize with new settings
}

void ADXL375_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin)
{
	ADXL_port = cs_port;
    ADXL_pin = cs_pin;
	ADXL_hspi = hspi;

	SPI_SwitchToMode3();

//	SPI_Write(hspi, cs_port, cs_pin, ADXL375_POWER_CTL, 0x00); // standby

    // Data format: Full resolution, ±200g (range bits = 00 for 200g)
//    SPI_Write(hspi, cs_port, cs_pin, ADXL375_DATA_FORMAT, 0x04);

    // Set bandwidth to 800 Hz (example)
    SPI_Write(hspi, cs_port, cs_pin, ADXL375_BW_RATE, 0x0F);

    // Measurement mode
    SPI_Write(hspi, cs_port, cs_pin, ADXL375_POWER_CTL, 0x08);

    SPI_RestoreMode();

    HAL_Delay(10);
}

void ADXL375_Read(Telemetry_t *telemetry)
{
	SPI_SwitchToMode3();

	uint8_t buffer[6];

    SPI_Read_Multi(ADXL_hspi, ADXL_port, ADXL_pin, ADXL375_DATAX0, buffer, 6);

    int16_t accel_x = (int16_t)(buffer[1] << 8 | buffer[0]);
    int16_t accel_y = (int16_t)(buffer[3] << 8 | buffer[2]);
    int16_t accel_z = (int16_t)(buffer[5] << 8 | buffer[4]);

    telemetry->adxl_accel_r = accel_x * 0.4805f;   // 0.049g * 9.80665
    telemetry->adxl_accel_p = -accel_y * 0.4805f;   // 0.049g * 9.80665
    telemetry->adxl_accel_y = -accel_z * 0.4805f;   // 0.049g * 9.80665

    SPI_RestoreMode();
}

uint8_t ADXL375_WhoAmI(void) {
	SPI_SwitchToMode3();

    uint8_t id = 0;
    SPI_Read(ADXL_hspi, ADXL_port, ADXL_pin, ADXL375_DEVID, &id, 1);

    SPI_RestoreMode();

    return id;
}

void LIS3MDLTR_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin) {
    LIS_port = cs_port;
    LIS_pin = cs_pin;
    LIS_hspi = hspi;

    // 1. Reset the device (CTRL_REG2: soft reset)
    SPI_Write(hspi, cs_port, cs_pin, LIS3MDLTR_CTRL_REG2, 0x0C); // soft reset + reboot
    HAL_Delay(50);

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

static float BMP388_compensate_temperature(uint32_t uncomp_temp, BMP388_CalibData *calib_data)
{
    float partial_data1;
    float partial_data2;

    partial_data1 = (float)(uncomp_temp - calib_data->par_t1);
    partial_data2 = (float)(partial_data1 * calib_data->par_t2);

    calib_data->t_lin = partial_data2 +
                        (partial_data1 * partial_data1) * calib_data->par_t3;

    return calib_data->t_lin;
}

static float BMP388_compensate_pressure(uint32_t uncomp_press, BMP388_CalibData *calib_data)
{
    float comp_press;
    float partial_data1;
    float partial_data2;
    float partial_data3;
    float partial_data4;
    float partial_out1;
    float partial_out2;

    partial_data1 = calib_data->par_p6 * calib_data->t_lin;
    partial_data2 = calib_data->par_p7 * (calib_data->t_lin * calib_data->t_lin);
    partial_data3 = calib_data->par_p8 * (calib_data->t_lin * calib_data->t_lin * calib_data->t_lin);
    partial_out1 = calib_data->par_p5 + partial_data1 + partial_data2 + partial_data3;

    partial_data1 = calib_data->par_p2 * calib_data->t_lin;
    partial_data2 = calib_data->par_p3 * (calib_data->t_lin * calib_data->t_lin);
    partial_data3 = calib_data->par_p4 * (calib_data->t_lin * calib_data->t_lin * calib_data->t_lin);
    partial_out2 = (float)uncomp_press *
                   (calib_data->par_p1 + partial_data1 + partial_data2 + partial_data3);

    partial_data1 = (float)uncomp_press * (float)uncomp_press;
    partial_data2 = calib_data->par_p9 + calib_data->par_p10 * calib_data->t_lin;
    partial_data3 = partial_data1 * partial_data2;

    partial_data4 = partial_data3 +
                    ((float)uncomp_press * (float)uncomp_press * (float)uncomp_press) *
                    calib_data->par_p11;

    comp_press = partial_out1 + partial_out2 + partial_data4;

    return comp_press;
}

void BMP388_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin)
{
    BMP388_port = cs_port;
    BMP388_pin = cs_pin;
    BMP388_hspi = hspi;
    HAL_Delay(20);

    SPI_Write(hspi, cs_port, cs_pin, BMP388_CMD, 0xB6);
    HAL_Delay(10);

    SPI_Write(hspi, cs_port, cs_pin, BMP388_PWR_CTRL, 0x33);

    SPI_Write(hspi, cs_port, cs_pin, BMP388_OSR, 0x03);

    SPI_Write(hspi, cs_port, cs_pin, BMP388_CONFIG, 0x02);

    HAL_Delay(10);

    //calibration value read
    uint8_t calib[21];
    SPI_Read_Multi(BMP388_hspi, BMP388_port, BMP388_pin, BMP388_CALIB_DATA, calib, 21);

    // Raw calibration values
    uint16_t T1 = (calib[1] << 8) | calib[0];
    uint16_t T2 = (calib[3] << 8) | calib[2];
    int8_t   T3 = calib[4];

    int16_t P1  = (calib[6] << 8) | calib[5];
    int16_t P2  = (calib[8] << 8) | calib[7];
    int8_t  P3  = calib[9];
    int8_t  P4  = calib[10];
    uint16_t P5 = (calib[12] << 8) | calib[11];
    uint16_t P6 = (calib[14] << 8) | calib[13];
    int8_t  P7  = calib[15];
    int8_t  P8  = calib[16];
    int16_t P9  = (calib[18] << 8) | calib[17];
    int8_t  P10 = calib[19];
    int8_t  P11 = calib[20];

    //convert based on datasheet
    bmp388_calib.par_t1 = T1 * 256.0f;
    bmp388_calib.par_t2 = T2 / 1073741824.0f;
    bmp388_calib.par_t3 = T3 / 281474976710656.0f;

    bmp388_calib.par_p1 = (P1 - 16384.0f) / 1048576.0f;
    bmp388_calib.par_p2 = (P2 - 16384.0f) / 536870912.0f;
    bmp388_calib.par_p3 = P3 / 4294967296.0f;
    bmp388_calib.par_p4 = P4 / 137438953472.0f;
    bmp388_calib.par_p5 = P5 * 8.0f;
    bmp388_calib.par_p6 = P6 / 64.0f;
    bmp388_calib.par_p7 = P7 / 256.0f;
    bmp388_calib.par_p8 = P8 / 32768.0f;
    bmp388_calib.par_p9 = P9 / 281474976710656.0f;
    bmp388_calib.par_p10 = P10 / 281474976710656.0f;
    bmp388_calib.par_p11 = P11 / 36893488147419103232.0f;
}

void BMP388_Read(Telemetry_t *telemetry)
{
    uint8_t buf[6];
    SPI_Read_Multi(BMP388_hspi, BMP388_port, BMP388_pin, BMP388_PRESS_DATA, buf, 6);

    int32_t raw_p = ((int32_t)buf[2] << 16) | ((int32_t)buf[1] << 8) | buf[0];
    int32_t raw_t = ((int32_t)buf[5] << 16) | ((int32_t)buf[4] << 8) | buf[3];

    //need to calibrate a lot
    telemetry->temperature2 = BMP388_compensate_temperature(raw_t, &bmp388_calib);
    telemetry->pressure2 = BMP388_compensate_pressure(raw_p, &bmp388_calib)/100.0f; // Pa → hPa

    telemetry->altitude2 = Calculate_Altitude(telemetry->pressure2);
}

uint8_t BMP388_WhoAmI(void)
{
    uint8_t id = 0;
    SPI_Read(BMP388_hspi, BMP388_port, BMP388_pin, BMP388_CHIP_ID, &id, 1);
    return id;
}
void calibrate_mag(Telemetry_t *telemetry){
	LIS3MDLTR_Read(telemetry);

	// Initialize max and min
	if (max_r == 0) max_r = telemetry->mag_r;
	if (max_p == 0) max_p = telemetry->mag_p;
	if (max_y == 0) max_y = telemetry->mag_y;
	if (min_r == 0) min_r = telemetry->mag_r;
	if (min_y == 0) min_y = telemetry->mag_y;
	if (min_p == 0) min_p = telemetry->mag_p;

	// Update max and min
	if (telemetry->mag_r > max_r) max_r = telemetry->mag_r;
	if (telemetry->mag_p > max_p) max_p = telemetry->mag_p;
	if (telemetry->mag_y > max_y) max_y = telemetry->mag_y;
	if (telemetry->mag_r < min_r) min_r = telemetry->mag_r;
	if (telemetry->mag_p < min_p) min_p = telemetry->mag_p;
	if (telemetry->mag_y < min_y) min_y = telemetry->mag_y;

	float offset_r = (max_r + min_r) / 2.0f;
	float offset_p = (max_p + min_p) / 2.0f;
	float offset_y = (max_y + min_y) / 2.0f;

	float radius_r = (max_r - min_r) / 2.0f;
	float radius_p = (max_p - min_p) / 2.0f;
	float radius_y = (max_y - min_y) / 2.0f;

	float avg_radius = (radius_r + radius_p + radius_y) / 3.0f;

	float scale_r = avg_radius / radius_r;
	float scale_p = avg_radius / radius_p;
	float scale_y = avg_radius / radius_y;

	write_mag("MAG_CALIB.csv", offset_r, offset_p, offset_y, scale_r, scale_p, scale_y);
}

void transform_accel_to_world(Telemetry_t *telemetry) {
  // Average IMUs (body frame)
  float ax = telemetry->lsm_accel_p;
  float ay = telemetry->lsm_accel_y;
  float az = telemetry->lsm_accel_r;

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
	HAL_GPIO_WritePin(IMU_CS_GPIO_Port, IMU_CS_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(IMU_2_CS_GPIO_Port, IMU_2_CS_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(Mag_CS_GPIO_Port, Mag_CS_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(Baro_CS_GPIO_Port, Baro_CS_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(Flash_CS_GPIO_Port, Flash_CS_Pin, GPIO_PIN_SET);
	HAL_GPIO_WritePin(SD_CS_GPIO_Port, SD_CS_Pin, GPIO_PIN_SET);
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

		sum_gr += telemetry->lsm_gyro_r + bias->lsm_gyro_r_bias;
		sum_gp += telemetry->lsm_gyro_p + bias->lsm_gyro_p_bias;
		sum_gy += telemetry->lsm_gyro_y + bias->lsm_gyro_y_bias;

		HAL_Delay(2);
	}

	// Divide by total samples to get average
	bias->lsm_gyro_r_bias = sum_gr / num_samples;
	bias->lsm_gyro_p_bias = sum_gp / num_samples;
	bias->lsm_gyro_y_bias = sum_gy / num_samples;
}

void Bias_Init(Bias_t *bias)
{
	bias->lsm_gyro_r_bias = -0.0784721822f;
	bias->lsm_gyro_p_bias = 0.0172601528f;
	bias->lsm_gyro_y_bias = 0.0110466592f;

	bias->lsm_accel_r_bias = 0.0f;
	bias->lsm_accel_p_bias = 0.0f;
	bias->lsm_accel_y_bias = 0.0f;

	bias->adxl_accel_r_bias = 0.0f;
	bias->adxl_accel_p_bias = 0.0f;
	bias->adxl_accel_y_bias = 0.0f;

	bias->bias_count = 0;

//    bias->adxl_accel_r_bias = 0.0f;
//    bias->adxl_accel_p_bias = 0.0f;
//    bias->adxl_accel_y_bias = 0.0f;
//
    bias->mag_r_bias = 25.816f; // 25.83, 25.305, 26.313
    bias->mag_p_bias = 2.196f; // 3.934, 1.547, 1.106
    bias->mag_y_bias = 20.629f; // 22.043, 21.798, 18.046

    bias->mag_r_scale = 1.0;
    bias->mag_p_scale = 1.0;
    bias->mag_y_scale = 1.0;

//    bias->bias_count = 0;
}

void Bias_Calculate(Bias_t *bias, Telemetry_t *t, int num_samples){

	Gyro_CalibrateBias(bias, t, num_samples);


	bias->bias_count = num_samples;

//    bias->bias_count += 1.0f;
//
//    float n = bias->bias_count;
//
//    bias->adxl_accel_r_bias += ((t->adxl_accel_r - 9.81) - bias->adxl_accel_r_bias) / n;
//    bias->adxl_accel_p_bias += (t->adxl_accel_p - bias->adxl_accel_p_bias) / n;
//    bias->adxl_accel_y_bias += (t->adxl_accel_y - bias->adxl_accel_y_bias) / n;
//
//    bias->lsm_accel_r_bias += ((t->lsm_accel_r - 9.81) - bias->lsm_accel_r_bias) / n;
//    bias->lsm_accel_p_bias += (t->lsm_accel_p - bias->lsm_accel_p_bias) / n;
//    bias->lsm_accel_y_bias += (t->lsm_accel_y - bias->lsm_accel_y_bias) / n;
}

void Apply_Bias(Bias_t *bias, Telemetry_t *t)
{

	t->lsm_gyro_r -= bias->lsm_gyro_r_bias;
	t->lsm_gyro_p -= bias->lsm_gyro_p_bias;
	t->lsm_gyro_y -= bias->lsm_gyro_y_bias;

	t->lsm_accel_r -= bias->lsm_accel_r_bias;
	t->lsm_accel_p -= bias->lsm_accel_p_bias;
	t->lsm_accel_y -= bias->lsm_accel_y_bias;

    t->adxl_accel_r -= bias->adxl_accel_r_bias;
    t->adxl_accel_p -= bias->adxl_accel_p_bias;
    t->adxl_accel_y -= bias->adxl_accel_y_bias;


    t->mag_r -= bias->mag_r_bias;
    t->mag_p -= bias->mag_p_bias;
    t->mag_y -= bias->mag_y_bias;

    t->mag_r *= bias->mag_r_scale;
    t->mag_p *= bias->mag_p_scale;
    t->mag_y *= bias->mag_y_scale;
}

// How to calibrate accel: run in debugger. Hold still.
// Place the board in 6 different orientations, letting gravity (1g) act on each axis positively and negatively.
// +X, -X, +Y, -Y, +Z, -Z facing down.
// For each position, record the average output of all three axes.
// For each axis, the offset is: Offset = (Value_+1g + Value_-1g) / 2
// The scale factor is: Scale = (Value_+1g - Value_-1g) / 2 (Theoretically, this should be 1g, but you can use it to correct minor gain errors).
float lsm_accel_r_bias_instance;
float lsm_accel_p_bias_instance;
float lsm_accel_y_bias_instance;
uint32_t prev_time_accel_bias;
void calibrate_accel_bias_stationary(Telemetry_t *telemetry)
{
	// Average IMUs (body frame)
	float ax = telemetry->lsm_accel_p;
	float ay = telemetry->lsm_accel_y;
	float az = telemetry->lsm_accel_r;

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
    lsm_accel_p_bias_instance = lsm_accel_p_bias_instance * alpha + bias_x * (alpha - 1);
    lsm_accel_y_bias_instance = lsm_accel_y_bias_instance * alpha + bias_y * (alpha - 1);
    lsm_accel_r_bias_instance = lsm_accel_r_bias_instance * alpha + bias_z * (alpha - 1);
}
