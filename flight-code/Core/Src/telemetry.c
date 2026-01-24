#include "telemetry.h"
#include "sd_card.h"
#include <stdio.h>

void init_data_file(SPI_HandleTypeDef *hspi){
	init_sd(hspi);
	write_sd(FLIGHT_DATA_FILE, "FLIGHT BEGIN");
	write_headers();
}

void log_data(char* state_str, Telemetry_t *t){
	char data_string[700];
	snprintf(data_string, sizeof(data_string),
	        "%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f",
	        state_str,
			t->pressure,
	        t->altitude,
	        t->startAlt,
	        t->temperature,
//	        t->angle_of_attack,
//	        t->velocity_r,
//	        t->velocity_p,
//	        t->velocity_y,
	        t->lsm_accel_r,
	        t->lsm_accel_p,
	        t->lsm_accel_y,
	        t->lsm_gyro_r,
	        t->lsm_gyro_p,
	        t->lsm_gyro_y
//	        t->icm_accel_r,
//	        t->icm_accel_p,
//	        t->icm_accel_y,
//	        t->icm_gyro_r,
//	        t->icm_gyro_p,
//	        t->icm_gyro_y,
//	        t->predicted_apogee,
//	        t->airbrake_deployment,
//	        t->mag_r,
//	        t->mag_p,
//	        t->mag_y
	    );

	write_sd(FLIGHT_DATA_FILE, data_string);
}

FRESULT write_headers(void)
{
    const char *header =
//        "state,pressure,altitude,startAlt,temperature,angle_of_attack,"
//        "velocity_r,velocity_p,velocity_y,"
//        "lsm_accel_r,lsm_accel_p,lsm_accel_y,"
//        "lsm_gyro_r,lsm_gyro_p,lsm_gyro_y,"
//        "icm_accel_r,icm_accel_p,icm_accel_y,"
//        "icm_gyro_r,icm_gyro_p,icm_gyro_y,"
//        "predicted_apogee,airbrake_deployment,"
//        "mag_r,mag_p,mag_y";
    	  "state,pressure,altitude,startAlt,temperature,"
          "lsm_accel_r,lsm_accel_p,lsm_accel_y,"
          "lsm_gyro_r,lsm_gyro_p,lsm_gyro_y";

    return write_sd(FLIGHT_DATA_FILE, header);
}
