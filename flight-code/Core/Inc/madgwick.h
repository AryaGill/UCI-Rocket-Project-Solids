//#pragma once
//#include "main.h"
//


#ifndef MADGWICK_H
#define MADGWICK_H

#include "telemetry.h"
#include <stdint.h>

void Madgwick_Init(Telemetry_t* telemetry, float beta);
void Madgwick_Update(Telemetry_t* telemetry);
void Madgwick_GetEuler(Telemetry_t* telemetry);

#endif // MADGWICK_H
