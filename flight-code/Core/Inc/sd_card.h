#pragma once

#include "stm32h7xx_hal.h"

void mount_sd();
void write_sd(char* file_name, char* data);
void read_sd(char* file_name, uint32_t line_number, char* RW_buffer, size_t buffer_size);
void delete_file(char* file_name);
void unmount_sd(char* file_name);
