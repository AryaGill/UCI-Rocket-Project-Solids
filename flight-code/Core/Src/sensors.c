#include "sensors.c"

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
	SPI_CS_HIGH(IMU_CS_GPIO_Port, IMU_CS_PIN);
	SPI_CS_HIGH(IMU_2_CS_GPIO_Port, IMU_2_CS_PIN);
	SPI_CS_HIGH(Baro_CS_GPIO_Port, Baro_CS_PIN);
	SPI_CS_HIGH(Mag_CS_GPIO_Port, Mag_CS_PIN);

	// Call init functions
	LSM6DSL_Init(hspi, IMU_CS_GPIO_Port, IMU_CS_PIN);
	ICM45686_Init(hspi, IMU_2_CS_GPIO_Port, IMU_2_CS_PIN);
	LPS22HH_Init(hspi, Baro_CS_GPIO_Port, Baro_CS_PIN);
	IIS2MDCTR_Init(hspi, Mag_CS_GPIO_Port, Mag_CS_PIN);
}

void read_sensors(SPI_HandleTypeDef *hspi, Telemetry_t *telemetry){
	LSM6DSL_Read(hspi, IMU_CS_GPIO_Port, IMU_CS_PIN, telemetry);
	ICM45686_Read(hspi, IMU_2_CS_GPIO_Port, IMU_2_CS_PIN, telemetry);
	LPS22HH_Read(hspi, Baro_CS_GPIO_Port, Baro_CS_PIN, telemetry);
	IIS2MDCTR_Read(hspi, Mag_CS_GPIO_Port, Mag_CS_PIN, telemetry);
}

void LSM6DSL_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin)
{
    SPI_Write(hspi, cs_port, cs_pin, LSM6DSL_CTRL3_C, 0x01); // reset
    HAL_Delay(10);

    SPI_Write(hspi, cs_port, cs_pin, LSM6DSL_CTRL1_XL, 0x60); // 416Hz, ±2g
    SPI_Write(hspi, cs_port, cs_pin, LSM6DSL_CTRL2_G,  0x60); // 416Hz, 245dps
}

void LSM6DSL_Read(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin, Telemetry_t *telemetry)
{
    uint8_t buf[12];
    SPI_Read(hspi, cs_port, cs_pin, LSM6DSL_OUTX_L_G, buf, 12);

    telemetry->gyro_r = buf[1]<<8 | buf[0]; // x
    telemetry->gyro_p = buf[3]<<8 | buf[2]; // y
    telemetry->gyro_y = buf[5]<<8 | buf[4]; // z
    telemetry->accel_r = buf[7]<<8 | buf[6]; // x
    telemetry->accel_p = buf[9]<<8 | buf[8]; // y
    telemetry->accel_y = buf[11]<<8 | buf[10]; // z
}

void ICM45686_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin)
{
    SPI_Write(hspi, cs_port, cs_pin, ICM_PWR_MGMT0, 0x0F); // accel+gyro ON
    HAL_Delay(5);

    SPI_Write(hspi, cs_port, cs_pin, ICM_GYRO_CFG0,  0x06); // 1kHz, ±500 dps
    SPI_Write(hspi, cs_port, cs_pin, ICM_ACCEL_CFG0, 0x06); // 1kHz, ±4g
}

void ICM45686_Read(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin, Telemetry_t *telemetry)
{
    uint8_t buf[12];
    SPI_Read(hspi, cs_port, cs_pin, ICM_DATA_START, buf, 12);


//    telemetry->ax = buf[1]<<8 | buf[0];
//    telemetry->ay = buf[3]<<8 | buf[2];
//    telemetry->az = buf[5]<<8 | buf[4];
//    telemetry->gx = buf[7]<<8 | buf[6];
//    telemetry->gy = buf[9]<<8 | buf[8];
//    telemetry->gz = buf[11]<<8 | buf[10];
}

void LPS22HH_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin)
{
    SPI_Write(hspi, cs_port, cs_pin, LPS22HH_CTRL_REG1, 0x50); // 75Hz, active
}

void LPS22HH_Read(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin, Telemetry_t *telemetry)
{
    uint8_t buf[5];
    SPI_Read(hspi, cs_port, cs_pin, LPS22HH_PRESS_OUT, buf, 5);

    telemetry->pressure = (int32_t)(buf[2]<<16 | buf[1]<<8 | buf[0]);

    telemetry->temperature = (int16_t)(buf[4]<<8 | buf[3]);
}

void IIS2MDCTR_Init(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin)
{
    SPI_Write(hspi, cs_port, cs_pin, IIS2M_CTRL_REG1, 0x9C); // 100Hz, continuous
}

void IIS2MDCTR_Read(SPI_HandleTypeDef *hspi, GPIO_TypeDef *cs_port, uint16_t cs_pin, Telemetry_t *telemetry)
{
    uint8_t buf[6];
    SPI_Read(hspi, cs_port, cs_pin, IIS2M_OUTX_L, buf, 6);

    telemetry->mag_r = buf[1]<<8 | buf[0]; // x
    telemetry->mag_p = buf[3]<<8 | buf[2]; // y
    telemetry->mag_y = buf[5]<<8 | buf[4]; // z
}
