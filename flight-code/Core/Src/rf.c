#include "rf.h"
#include <stdio.h>
#include <string.h>

static UART_HandleTypeDef *rf_huart;

void RF_Init(UART_HandleTypeDef *huart) {
    rf_huart = huart;
}

void RF_Transmit(Telemetry_t *telemetry) {
    char tx_buffer[128];

    // Format: "P:%.2f,T:%.2f,A:%.2f\r\n"
    // P = Pressure (hPa), T = Temperature (°C), A = Altitude (m)
    int len = snprintf(tx_buffer, sizeof(tx_buffer),
                      "P:%.2f,T:%.2f,A:%.2f\r\n",
                      telemetry->pressure,
                      telemetry->temperature,
                      telemetry->altitude);

    // Transmit over UART4 (Serial 4)
    HAL_UART_Transmit(rf_huart, (uint8_t*)tx_buffer, len, HAL_MAX_DELAY);
}
