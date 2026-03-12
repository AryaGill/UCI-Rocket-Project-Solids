#include "sd_card.h"
#include "main.h"
#include "telemetry.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// SPI Helper Functions
static inline void SPI_CS_LOW(GPIO_TypeDef *port, uint16_t pin)
{
    HAL_GPIO_WritePin(port, pin, GPIO_PIN_RESET);
}

static inline void SPI_CS_HIGH(GPIO_TypeDef *port, uint16_t pin)
{
    HAL_GPIO_WritePin(port, pin, GPIO_PIN_SET);
}

void SD_SendDummyClocks(SPI_HandleTypeDef *hspi,
                        GPIO_TypeDef *port, uint16_t pin)
{
    uint8_t dummy = 0xFF;

    SPI_CS_HIGH(port, pin);

    for (int i = 0; i < 10; i++) { // 10 bytes = 80 clocks
        HAL_SPI_Transmit(hspi, &dummy, 1, HAL_MAX_DELAY);
    }
}

void init_sd(SPI_HandleTypeDef *hspi){
	SD_SendDummyClocks(hspi, SD_CS_GPIO_Port, SD_CS_Pin);

	SPI_CS_LOW(SD_CS_GPIO_Port, SD_CS_Pin);
	f_mount(&USERFatFS, USERPath, 1);
	SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
}

FRESULT open_file(File_t *f){
	FRESULT res;

	SPI_CS_LOW(SD_CS_GPIO_Port, SD_CS_Pin);

	// Open file (create if it doesn't exist)
	res = f_open(&f->file, f->file_name, FA_OPEN_ALWAYS | FA_WRITE);
	if (res != FR_OK){
		SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
		return res;
	}

	// Move write pointer to end of file
	res = f_lseek(&f->file, f_size(&f->file));
	if (res != FR_OK) {
		f_close(&f->file);
		SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
		return res;
	}

	SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);

	f->file_open = 1;
	return FR_OK;
}

FRESULT write_sd(File_t *f, const char *line)
{
	FRESULT res;

	if (!f->file_open){
		res = open_file(f);
		if (res != FR_OK){
			return res;
		}
	}

    UINT bytes_written;
    SPI_CS_LOW(SD_CS_GPIO_Port, SD_CS_Pin);

    // Write the line
    res = f_write(&f->file, line, strlen(line), &bytes_written);
    if (res != FR_OK) {
        f_close(&f->file);
        SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
        return res;
    }

    // newline
    const char newline[] = "\r\n";
    f_write(&f->file, newline, 2, &bytes_written);

    SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
    return FR_OK;
}

FRESULT flush_file(File_t *f)
{
    if (!f->file_open)
        return FR_INVALID_OBJECT;

    return f_sync(&f->file);
}

//FRESULT read_sd_line(const char *filename, char *buffer, UINT buffer_size)
//{
//    FIL file;
//    FRESULT res;
//
//    // Open file for reading
//    res = f_open(&file, filename, FA_READ);
//    if (res != FR_OK)
//        return res;
//
//    // Read one line
//    if (f_gets(buffer, buffer_size, &file) == NULL) {
//        f_close(&file);
//        return FR_DISK_ERR;   // or FR_DISK_ERR if you prefer
//    }
//
//    f_close(&file);
//    return FR_OK;
//}

FRESULT write_sd_state(const char *filename, FlightState_t state, float start_alt){
	FIL file;
	UINT bytes_written;
	FRESULT res;

	SPI_CS_LOW(SD_CS_GPIO_Port, SD_CS_Pin);

	// Open file: create new or overwrite existing
	res = f_open(&file, filename, FA_WRITE | FA_CREATE_ALWAYS);
	if (res != FR_OK){
		SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
		return res;
	}

	// Write state
	char line[20];
	state_to_string_num(state, line);
	res = f_write(&file, line, strlen(line), &bytes_written);
	if (res != FR_OK) {
		f_close(&file);
		SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
		return res;
	}

	// newline
	const char newline[] = "\r\n";
	f_write(&file, newline, 2, &bytes_written);

	// Write start alt
	snprintf(line, sizeof(line), "%.6f", start_alt);
	res = f_write(&file, line, strlen(line), &bytes_written);
	if (res != FR_OK) {
		f_close(&file);
		SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
		return res;
	}

	// newline
	f_write(&file, newline, 2, &bytes_written);

	f_close(&file);

	SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
	return FR_OK;
}

void strip_newline(char *s)
{
    size_t len = strlen(s);
    while (len > 0 &&
          (s[len-1] == '\n' || s[len-1] == '\r'))
    {
        s[--len] = '\0';
    }
}

FRESULT read_sd_state(const char *filename, FlightState_t *state, float *start_alt){
	FIL file;
	FRESULT res;

	char buf[20];

	SPI_CS_LOW(SD_CS_GPIO_Port, SD_CS_Pin);

	// Open file for reading
	res = f_open(&file, filename, FA_READ);
	if (res != FR_OK){
		SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
		return res;
	}

	// Read state line
	if (f_gets(buf, sizeof(buf), &file) == NULL) {
		f_close(&file);
		SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
		return FR_DISK_ERR;   // or FR_DISK_ERR if you prefer
	}
	strip_newline(buf);
	int state_int = atoi(buf);
	if (state_int < (int)DISARMED || state_int > (int)LANDED){
		*state = DISARMED;
	}
	else{
		*state = (FlightState_t)atoi(buf);
	}

	// Read start alt line
	if (f_gets(buf, sizeof(buf), &file) == NULL) {
		f_close(&file);
		SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
		return FR_DISK_ERR;   // or FR_DISK_ERR if you prefer
	}
	strip_newline(buf);
	*start_alt = strtof(buf, NULL);

	f_close(&file);

	SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
	return FR_OK;
}

uint8_t sd_file_exists(const char *filename)
{
    FILINFO fno;
    FRESULT res;

    SPI_CS_LOW(SD_CS_GPIO_Port, SD_CS_Pin);

    res = f_stat(filename, &fno);

    SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);

    if (res == FR_OK)
        return 1;   // file exists
    else
        return 0;   // file does not exist or error
}

FRESULT sd_delete_file(const char *filename) {
	SPI_CS_LOW(SD_CS_GPIO_Port, SD_CS_Pin);
	FRESULT res = f_unlink(filename);
	SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
    return res;
}

FRESULT write_mag(const char *filename,
                  float mag_r_bias, float mag_r_scale,
                  float mag_p_bias, float mag_p_scale,
                  float mag_y_bias, float mag_y_scale)
{
    FIL file;
    UINT bytes_written;
    char line[120];
    FRESULT res;

    SPI_CS_LOW(SD_CS_GPIO_Port, SD_CS_Pin);

    // Open file and overwrite if it exists
    res = f_open(&file, filename, FA_CREATE_ALWAYS | FA_WRITE);
    if (res != FR_OK) {
        SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
        return res;
    }

    snprintf(line, sizeof(line),
        "mag_r_bias: %.6f\n"
        "mag_r_scale: %.6f\n"
        "mag_p_bias: %.6f\n"
        "mag_p_scale: %.6f\n"
        "mag_y_bias: %.6f\n"
        "mag_y_scale: %.6f\n",
        mag_r_bias, mag_r_scale,
        mag_p_bias, mag_p_scale,
        mag_y_bias, mag_y_scale);

    res = f_write(&file, line, strlen(line), &bytes_written);
    if (res != FR_OK) {
        f_close(&file);
        SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
        return res;
    }

    f_close(&file);

    SPI_CS_HIGH(SD_CS_GPIO_Port, SD_CS_Pin);
    return FR_OK;
}
