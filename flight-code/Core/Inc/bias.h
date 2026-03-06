#pragma once

#ifndef BIAS_H
#define BIAS_H

#include "main.h"

void Bias_Init(Bias_t *bias);

void Bias_Calculate(Bias_t *bias, Telemetry_t *telemetry);

void Apply_Bias(Bias_t *bias, Telemetry_t *telemetry);

#endif
