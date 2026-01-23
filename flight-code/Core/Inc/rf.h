#ifndef RF_H
#define RF_H

#include "main.h"
#include "telemetry.h"

void RF_Init(UART_HandleTypeDef *huart);
void RF_Transmit(Telemetry_t *telemetry);

#endif // RF_H
