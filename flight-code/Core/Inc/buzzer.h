#pragma once
#include "main.h"

#define BUZZER_CHANNEL TIM_CHANNEL_1

#define BUZZER_ON 1 // Set to 1 for flight. Set to 0 for testing

// Turn on and set buzzer freq (Hz)
void buzzer_set_frequency(uint32_t freq);
