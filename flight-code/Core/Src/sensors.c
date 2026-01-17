#include "sensors.h"

GPIO_TypeDef *LPS22HH_port;
uint16_t LPS22HH_pin;
SPI_HandleTypeDef *LPS22HH_hspi;

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
    // Set CS pin high (inactive)
    SPI_CS_HIGH(Baro_CS_GPIO_Port, Baro_CS_Pin);
    HAL_Delay(10);

    // Initialize LPS22HHTR
    LPS22HH_Init(hspi, Baro_CS_GPIO_Port, Baro_CS_Pin);
}

// Sensor Reading
void read_sensors(Telemetry_t *telemetry)
{
    LPS22HH_Read(telemetry);
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
