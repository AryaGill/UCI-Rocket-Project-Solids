#include "telemetry.h"
#include "sd_card.h"
#include <stdio.h>

void init_data_file(SPI_HandleTypeDef *hspi){
	init_sd(hspi);
	write_sd(FLIGHT_DATA_FILE, "FLIGHT BEGIN");
	write_headers();
}

void log_data(FlightState_t flight_state, Telemetry_t *t){
	char data_string[700];
	char state_str[3];
	state_to_string_num(flight_state, state_str);
	snprintf(data_string, sizeof(data_string),
	        "%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%u,%u",
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
	        t->lsm_gyro_y,
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
//	        t->mag_y,
			t->q0,
			t->q1,
			t->q2,
			t->q3,
			t->accel_world_x,
			t->accel_world_y,
			t->accel_world_z,
			t->alt_fused,
			t->cam1_on,
			t->cam2_on
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

void state_to_string_num(FlightState_t state, char* str) {
	switch(state) {
		case DISARMED:
			strcpy(str, "0");
			break;
    	case LAUNCH_PAD:
    		strcpy(str, "1");
    		break;
    	case MOTOR_BURN:
    		strcpy(str, "2");
    		break;
    	case GLIDING_ASCENT:
    		strcpy(str, "3");
    		break;
    	case DROGUE_PRIMARY_DEPLOYING:
    		strcpy(str, "4");
    		break;
    	case DROGUE_PRIMARY_DEPLOYED:
    		strcpy(str, "5");
    		break;
    	case DROGUE_SECONDARY_DEPLOYING:
    		strcpy(str, "6");
    		break;
    	case DROGUE_SECONDARY_DEPLOYED:
    		strcpy(str, "7");
    		break;
    	case MAIN_PRIMARY_DEPLOYING:
    		strcpy(str, "8");
    		break;
    	case MAIN_PRIMARY_DEPLOYED:
    		strcpy(str, "9");
    		break;
    	case MAIN_SECONDARY_DEPLOYING:
    		strcpy(str, "10");
    		break;
    	case MAIN_SECONDARY_DEPLOYED:
    		strcpy(str, "11");
    		break;
    	case LANDED:
    		strcpy(str, "12");
    		break;
    	default:
    		strcpy(str, "-1");
    		break;
	}
}

void state_to_string_name(FlightState_t state, char* str) {
	switch(state) {
		case DISARMED:
			strcpy(str, "DISARMED");
			break;
    	case LAUNCH_PAD:
    		strcpy(str, "LAUNCH_PAD");
    		break;
    	case MOTOR_BURN:
    		strcpy(str, "MOTOR_BURN");
    		break;
    	case GLIDING_ASCENT:
    		strcpy(str, "GLIDING_ASCENT");
    		break;
    	case DROGUE_PRIMARY_DEPLOYING:
    		strcpy(str, "DROGUE_PRIMARY_DEPLOYING");
    		break;
    	case DROGUE_PRIMARY_DEPLOYED:
    		strcpy(str, "DROGUE_PRIMARY_DEPLOYED");
    		break;
    	case DROGUE_SECONDARY_DEPLOYING:
    		strcpy(str, "DROGUE_SECONDARY_DEPLOYING");
    		break;
    	case DROGUE_SECONDARY_DEPLOYED:
    		strcpy(str, "DROGUE_SECONDARY_DEPLOYED");
    		break;
    	case MAIN_PRIMARY_DEPLOYING:
    		strcpy(str, "MAIN_PRIMARY_DEPLOYING");
    		break;
    	case MAIN_PRIMARY_DEPLOYED:
    		strcpy(str, "MAIN_PRIMARY_DEPLOYED");
    		break;
    	case MAIN_SECONDARY_DEPLOYING:
    		strcpy(str, "MAIN_SECONDARY_DEPLOYING");
    		break;
    	case MAIN_SECONDARY_DEPLOYED:
    		strcpy(str, "MAIN_SECONDARY_DEPLOYED");
    		break;
    	case LANDED:
    		strcpy(str, "LANDED");
    		break;
    	default:
    		strcpy(str, "UNKNOWN");
    		break;
	}
}
