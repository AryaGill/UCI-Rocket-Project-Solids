#include "sensors.h"

#define LSM6DSL_WHO_AM_I 0x0F //temp remove later

GPIO_TypeDef *LSM6DSL_port;
uint16_t LSM6DSL_pin;
SPI_HandleTypeDef *LSM6DSL_hspi;

GPIO_TypeDef *ICM45686_port;
uint16_t ICM45686_pin;
SPI_HandleTypeDef *ICM45686_hspi;

GPIO_TypeDef *LPS22HH_port;
uint16_t LPS22HH_pin;
SPI_HandleTypeDef *LPS22HH_hspi;

GPIO_TypeDef *IIS2MDCTR_port;
uint16_t IIS2MDCTR_pin;
SPI_HandleTypeDef *IIS2MDCTR_hspi;

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
    reg |= 0x80; // read
    SPI_CS_LOW(port, pin);
    HAL_SPI_Transmit(hspi, &reg, 1, HAL_MAX_DELAY);
    HAL_SPI_Receive(hspi, buf, len, HAL_MAX_DELAY);
    SPI_CS_HIGH(port, pin);
}

static void SPI_Write(SPI_HandleTypeDef *hspi, GPIO_TypeDef *port, uint16_t pin,
                      uint8_t reg, uint8_t val)
{
    uint8_t tx[2] = { reg & 0x7F, val };
    SPI_CS_LOW(port, pin);
    HAL_SPI_Transmit(hspi, tx, 2, HAL_MAX_DELAY);
    SPI_CS_HIGH(port, pin);
}

void init_sensors(SPI_HandleTypeDef *hspi){
	// Set all CS pins high
	SPI_CS_HIGH(IMU_CS_GPIO_Port, IMU_CS_Pin);
	SPI_CS_HIGH(IMU_2_CS_GPIO_Port, IMU_2_CS_Pin);
	SPI_CS_HIGH(Baro_CS_GPIO_Port, Baro_CS_Pin);
	SPI_CS_HIGH(Mag_CS_GPIO_Port, Mag_CS_Pin);

	// Call init functions
	LSM6DSL_Init(hspi, IMU_CS_GPIO_Port, IMU_CS_Pin);
	ICM45686_Init(hspi, IMU_2_CS_GPIO_Port, IMU_2_CS_Pin);
	LPS22HH_Init(hspi, Baro_CS_GPIO_Port, Baro_CS_Pin);
	IIS2MDCTR_Init(hspi, Mag_CS_GPIO_Port, Mag_CS_Pin);
}

void read_sensors(Telemetry_t *telemetry){
	LSM6DSL_Read(telemetry);
	ICM45686_Read(telemetry);
	LPS22HH_Read(telemetry);
	IIS2MDCTR_Read(telemetry);
}

void LSM6DSL_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin)
{
	LSM6DSL_port = cs_port;
	LSM6DSL_pin = cs_pin;
	LSM6DSL_hspi = hspi;

    SPI_Write(hspi, cs_port, cs_pin, LSM6DSL_CTRL3_C, 0x01); // reset
    HAL_Delay(10);

    SPI_Write(hspi, cs_port, cs_pin, LSM6DSL_CTRL1_XL, 0x60); // 416Hz, ±2g
    SPI_Write(hspi, cs_port, cs_pin, LSM6DSL_CTRL2_G,  0x60); // 416Hz, 245dps
}

void LSM6DSL_Read(Telemetry_t *telemetry)
{
    uint8_t buf[12];
    int16_t gx, gy, gz;
    int16_t ax, ay, az;

    SPI_Read(LSM6DSL_hspi, LSM6DSL_port, LSM6DSL_pin, LSM6DSL_OUTX_L_G, buf, 12);

    /* Raw signed values */
	gx = (int16_t)(buf[1]  << 8 | buf[0]);
	gy = (int16_t)(buf[3]  << 8 | buf[2]);
	gz = (int16_t)(buf[5]  << 8 | buf[4]);
	ax = (int16_t)(buf[7]  << 8 | buf[6]);
	ay = (int16_t)(buf[9]  << 8 | buf[8]);
	az = (int16_t)(buf[11] << 8 | buf[10]);

	/* Scale to metric units */
	telemetry->gyro_r  = gx * 0.000152716f;   // rad/s
	telemetry->gyro_p  = gy * 0.000152716f;   // rad/s
	telemetry->gyro_y  = gz * 0.000152716f;   // rad/s

	telemetry->accel_r = ax * 0.00059855f;    // m/s²
	telemetry->accel_p = ay * 0.00059855f;    // m/s²
	telemetry->accel_y = az * 0.00059855f;    // m/s²
}

void ICM45686_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin)
{
	ICM45686_port = cs_port;
	ICM45686_pin = cs_pin;
	ICM45686_hspi = hspi;

	SPI_Write(hspi, cs_port, cs_pin, ICM_PWR_MGMT0, 0x0F); // accel+gyro ON
    HAL_Delay(5);

    SPI_Write(hspi, cs_port, cs_pin, ICM_GYRO_CFG0,  0x06); // 1kHz, ±500 dps
    SPI_Write(hspi, cs_port, cs_pin, ICM_ACCEL_CFG0, 0x06); // 1kHz, ±4g
}

void ICM45686_Read(Telemetry_t *telemetry)
{
	uint8_t buf[12];
	int16_t ax, ay, az;
	int16_t gx, gy, gz;

	SPI_Read(ICM45686_hspi, ICM45686_port, ICM45686_pin, ICM_DATA_START, buf, 12);

	/* Raw signed values */
	ax = (int16_t)(buf[1]  << 8 | buf[0]);
	ay = (int16_t)(buf[3]  << 8 | buf[2]);
	az = (int16_t)(buf[5]  << 8 | buf[4]);
	gx = (int16_t)(buf[7]  << 8 | buf[6]);
	gy = (int16_t)(buf[9]  << 8 | buf[8]);
	gz = (int16_t)(buf[11] << 8 | buf[10]);

	/* Convert to metric units */
	telemetry->accel_r = ax * 0.001197f;     // m/s²
	telemetry->accel_p = ay * 0.001197f;     // m/s²
	telemetry->accel_y = az * 0.001197f;     // m/s²

	telemetry->gyro_r = gx * 0.00026646f;   // rad/s
	telemetry->gyro_p = gy * 0.00026646f;   // rad/s
	telemetry->gyro_y = gz * 0.00026646f;   // rad/s
}

void LPS22HH_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin)
{
	LPS22HH_port = cs_port;
	LPS22HH_pin = cs_pin;
	LPS22HH_hspi = hspi;

	SPI_Write(hspi, cs_port, cs_pin, LPS22HH_CTRL_REG1, 0x50); // 75Hz, active
}

void LPS22HH_Read(Telemetry_t *telemetry)
{
	uint8_t buf[5];
	int32_t raw_p;
	int16_t raw_t;

	SPI_Read(LPS22HH_hspi, LPS22HH_port, LPS22HH_pin, LPS22HH_PRESS_OUT, buf, 5);

	/* Raw values */
	raw_p = (int32_t)(buf[2] << 16 | buf[1] << 8 | buf[0]);
	raw_t = (int16_t)(buf[4] << 8 | buf[3]);

	/* Convert to metric */
	telemetry->pressure = raw_p / 4096.0f;   // hPa
	telemetry->temperature = raw_t / 100.0f; // °C
}

// Set SDIO direction
static inline void MAG_SDIO_OUT(GPIO_TypeDef *port, uint16_t pin) {
    GPIO_InitTypeDef GPIO_InitStruct = {
        .Pin = pin,
        .Mode = GPIO_MODE_OUTPUT_PP,
        .Pull = GPIO_NOPULL,
        .Speed = GPIO_SPEED_FREQ_HIGH
    };
    HAL_GPIO_Init(port, &GPIO_InitStruct);
}

static inline void MAG_SDIO_IN(GPIO_TypeDef *port, uint16_t pin) {
    GPIO_InitTypeDef GPIO_InitStruct = {
        .Pin = pin,
        .Mode = GPIO_MODE_INPUT,
        .Pull = GPIO_NOPULL
    };
    HAL_GPIO_Init(port, &GPIO_InitStruct);
}

// Pulse SCLK once
static inline void SPI_SCLK_PULSE(SPI_HandleTypeDef *hspi) {
    uint8_t dummy = 0xFF;
    HAL_SPI_Transmit(hspi, &dummy, 1, HAL_MAX_DELAY); // Only clock edges needed
}

static void MAG_WriteReg(SPI_HandleTypeDef *hspi, GPIO_TypeDef *port, uint16_t pin, uint8_t reg, uint8_t val) {
    SPI_CS_LOW(port, pin);
    MAG_SDIO_OUT(port, pin);

    // MSB first, set write bit = 0
    for (int i = 7; i >= 0; i--) {
        HAL_GPIO_WritePin(port, pin,
                          (reg >> i) & 1 ? GPIO_PIN_SET : GPIO_PIN_RESET);
        SPI_SCLK_PULSE(hspi);
    }

    // Send value byte
    for (int i = 7; i >= 0; i--) {
        HAL_GPIO_WritePin(port, pin,
                          (val >> i) & 1 ? GPIO_PIN_SET : GPIO_PIN_RESET);
        SPI_SCLK_PULSE(hspi);
    }

    SPI_CS_HIGH(port, pin);
}

static uint8_t MAG_ReadReg(SPI_HandleTypeDef *hspi, GPIO_TypeDef *port, uint16_t pin, uint8_t reg) {
    uint8_t val = 0;
    SPI_CS_LOW(port, pin);
    MAG_SDIO_OUT(port, pin);

    // Set read bit = 1
    reg |= 0x80;

    // Send register address
    for (int i = 7; i >= 0; i--) {
        HAL_GPIO_WritePin(port, pin,
                          (reg >> i) & 1 ? GPIO_PIN_SET : GPIO_PIN_RESET);
        SPI_SCLK_PULSE(hspi);
    }

    // Switch to input for reading
    MAG_SDIO_IN(port, pin);

    // Read value
    for (int i = 7; i >= 0; i--) {
    	SPI_SCLK_PULSE(hspi);
        if (HAL_GPIO_ReadPin(port, pin))
            val |= (1 << i);
    }

    SPI_CS_HIGH(port, pin);
    return val;
}


void IIS2MDCTR_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin) {
	IIS2MDCTR_port = cs_port;
	IIS2MDCTR_pin = cs_pin;
	IIS2MDCTR_hspi = hspi;

	// CTRL_REG1 = 0x9C → continuous mode, 100Hz, SIM=1 (3-wire)
    MAG_WriteReg(hspi, cs_port, cs_pin, 0x20, 0x9C);
}

void IIS2MDCTR_Read(Telemetry_t *telemetry) {
    uint8_t buf[6];
    int16_t mx, my, mz;

    // Read 6 consecutive registers: OUTX_L, OUTX_H, OUTY_L, ...
    for (uint8_t i = 0; i < 6; i++) {
        buf[i] = MAG_ReadReg(IIS2MDCTR_hspi, IIS2MDCTR_port, IIS2MDCTR_pin, 0x28 + i);
    }

    // Combine bytes
    mx = (int16_t)(buf[1] << 8 | buf[0]);
    my = (int16_t)(buf[3] << 8 | buf[2]);
    mz = (int16_t)(buf[5] << 8 | buf[4]);

    // Scale to µT (1.5 mGauss/LSB → 0.15 µT)
    telemetry->mag_r = mx * 0.15f;
    telemetry->mag_p = my * 0.15f;
    telemetry->mag_y = mz * 0.15f;
}

uint8_t LSM6DSL_WhoAmI(void)
{
    uint8_t id = 0;
    SPI_Read(LSM6DSL_hspi,
             LSM6DSL_port,
             LSM6DSL_pin,
             LSM6DSL_WHO_AM_I,
             &id,
             1);
    return id;
}
uint8_t LSM6DSL_ReadReg(uint8_t reg)
{
    uint8_t val = 0;
    SPI_Read(LSM6DSL_hspi,
             LSM6DSL_port,
             LSM6DSL_pin,
             reg,
             &val,
             1);
    return val;
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

void LPS22HH_TestRead(float *pressure_hpa)
{
    uint8_t buf[3];
    int32_t raw_p;

    // Pressure output registers start at 0x28 (XL, L, H)
    SPI_Read(LPS22HH_hspi,
             LPS22HH_port,
             LPS22HH_pin,
             LPS22HH_PRESS_OUT,
             buf,
             3);

    raw_p = (int32_t)(buf[2] << 16 | buf[1] << 8 | buf[0]);

    // Datasheet: pressure = raw / 4096 hPa
    *pressure_hpa = raw_p / 4096.0f;
}
void ICM45686_TestRead(int16_t *ax, int16_t *ay, int16_t *az)
{
    uint8_t buf[6];

    // Read accel X/Y/Z (first 6 bytes of data block)
    SPI_Read(ICM45686_hspi,
             ICM45686_port,
             ICM45686_pin,
             ICM_DATA_START,
             buf,
             6);

    *ax = (int16_t)(buf[1] << 8 | buf[0]);
    *ay = (int16_t)(buf[3] << 8 | buf[2]);
    *az = (int16_t)(buf[5] << 8 | buf[4]);
}
void IIS2MDCTR_TestRead(int16_t *mx, int16_t *my, int16_t *mz)
{
    uint8_t buf[6];

    for (uint8_t i = 0; i < 6; i++) {
        buf[i] = MAG_ReadReg(IIS2MDCTR_hspi,
                             IIS2MDCTR_port,
                             IIS2MDCTR_pin,
                             IIS2M_OUTX_L + i);
    }

    *mx = (int16_t)(buf[1] << 8 | buf[0]);
    *my = (int16_t)(buf[3] << 8 | buf[2]);
    *mz = (int16_t)(buf[5] << 8 | buf[4]);
}


void LSM6DSL_WriteReg(uint8_t reg, uint8_t val)
{
    SPI_Write(LSM6DSL_hspi,
              LSM6DSL_port,
              LSM6DSL_pin,
              reg,
              val);
}

