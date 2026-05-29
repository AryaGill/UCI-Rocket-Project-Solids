#pragma once

#include "stm32h7xx_hal.h"
#include "fatfs.h"

typedef struct {
	FIL file;
	char file_name[20];
	uint8_t file_open;
} File_t;

void init_sd(SPI_HandleTypeDef *hspi);
FRESULT open_file(File_t *f);
FRESULT write_sd(File_t *f, const char *line);
FRESULT flush_file(File_t *f);
FRESULT write_sd_state(const char *filename, FlightState_t state, float start_alt);
FRESULT read_sd_state(const char *filename, FlightState_t *state, float *start_alt);
uint8_t sd_file_exists(const char *filename);
FRESULT sd_delete_file(const char *filename);
FRESULT sd_remove_file(File_t *f);

FRESULT write_mag(const char *filename, float mag_r_bias, float mag_r_scale, float mag_p_bias, float mag_p_scale, float mag_y_bias, float mag_y_scale);
