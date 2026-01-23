#include "rf.h"
#include <stdio.h>
#include <string.h>

static UART_HandleTypeDef *rf_huart;

void RF_Init(UART_HandleTypeDef *huart) {
    rf_huart = huart;
}

void RF_Transmit(Telemetry_t *telemetry) {
    char tx_buffer[128];

    // Pressure (hPa), Temperature (°C), Altitude (m)
    int len = snprintf(tx_buffer, sizeof(tx_buffer),
                      "%.2f,%.2f,%.2f\r\n",
                      telemetry->pressure,
                      telemetry->temperature,
                      telemetry->altitude);

    // Transmit over UART4 (Serial 4)
    HAL_UART_Transmit(rf_huart, (uint8_t*)tx_buffer, len, HAL_MAX_DELAY);
}
