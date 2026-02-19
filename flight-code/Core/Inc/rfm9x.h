#pragma once

#include "main.h"

/* Public API */

void RFM9X_Init(SPI_HandleTypeDef *hspi_p, GPIO_TypeDef *cs_port_p, uint16_t cs_pin_p, GPIO_TypeDef *reset_port_p, uint16_t reset_pin_p, GPIO_TypeDef *en_port_p, uint16_t en_pin_p);
void RFM9X_Reset();

void RFM9X_SetFrequency(uint32_t freq_hz);
void RFM9X_SetTxPower(uint8_t power);

void RFM9X_Send(uint8_t *data, uint8_t len);
void RFM9X_Poll();

uint8_t RFM9X_IsTxBusy();

uint8_t RFM9X_Receive(uint8_t *buf, uint8_t max_len);
