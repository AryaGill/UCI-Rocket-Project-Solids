#include "telemetry.h"
#include "sd_card.h"
#include <stdio.h>
#include <string.h>

File_t data_file = {0};

void init_data_file(SPI_HandleTypeDef *hspi){
	init_sd(hspi);
	strcpy(data_file.file_name, FLIGHT_DATA_FILE);
	open_file(&data_file);
	write_sd(&data_file, "FLIGHT BEGIN");
	write_headers();
	flush_file(&data_file);
}

void get_rf_msg(FlightState_t flight_state, Telemetry_t *t, char* msg, size_t msg_size){
	char state_str[3];
	state_to_string_num(flight_state, state_str);
	snprintf(msg, msg_size,
			"%lu,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%u,%u,%u,%lu,%lu,%lu,%lu,%s,%s\r\n",
			  t->time,
			  t->temperature,
			  t->pressure,
			  t->altitude - t->startAlt,
			  t->pressure, // todo: change back to t->alt_fused
			  t->lsm_gyro_r,
			  t->lsm_gyro_p,
			  t->lsm_gyro_y,
			  t->lsm_accel_r,
			  t->lsm_accel_p,
			  t->lsm_accel_y,
			  t->accel_world_x,
			  t->accel_world_y,
			  t->accel_world_z,
			  t->mag_r,
			  t->mag_p,
			  t->mag_y,
			  t->roll,
			  t->pitch,
			  t->yaw,
			  t->velocity_world_x,
			  t->velocity_world_y,
			  t->velocity_world_z,
			  t->q0,
			  t->q1,
			  t->q2,
			  t->q3,
			  t->predicted_apogee,
			  t->airbrake_deployment,
			  t->cam1_on,
			  t->cam2_on,
			  t->main_p_ematch_voltage,
			  t->main_s_ematch_voltage,
			  t->drogue_p_ematch_voltage,
			  t->drogue_s_ematch_voltage,
			  state_str,
			  t->cmd_echo);
}

void log_data(FlightState_t flight_state, Telemetry_t *t){
	char data_string[1000];
	char state_str[3];
	state_to_string_num(flight_state, state_str);
	snprintf(data_string, sizeof(data_string),
	        "%lu,%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%u,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%u,%u,%lu,%lu,%lu,%lu,%.3f,%.3f,%.3f,%lu,%lu,%lu,%lu,%lu,%i,%i",
	        t->time,
			state_str,
			t->pressure,
	        t->altitude,
	        t->startAlt,
	        t->temperature,
			t->velocity_world_x,
			t->velocity_world_y,
			t->velocity_world_z,
			t->baro_vz,
	        t->lsm_accel_r,
	        t->lsm_accel_p,
	        t->lsm_accel_y,
	        t->lsm_gyro_r,
	        t->lsm_gyro_p,
	        t->lsm_gyro_y,
	        t->adxl_accel_r,
	        t->adxl_accel_p,
	        t->adxl_accel_y,
	        t->predicted_apogee,
	        t->airbrake_deployment,
	        t->mag_r,
	        t->mag_p,
	        t->mag_y,
			t->q0,
			t->q1,
			t->q2,
			t->q3,
			t->accel_world_x,
			t->accel_world_y,
			t->accel_world_z,
			t->alt_fused,
			t->cam1_on,
			t->cam2_on,
			t->main_p_ematch_voltage,
			t->main_s_ematch_voltage,
			t->drogue_p_ematch_voltage,
			t->drogue_s_ematch_voltage,
			t->roll,
			t->pitch,
			t->yaw,
			t->t_burnout,
			t->t_apogee,
			t->t_drogue,
			t->t_main,
			t->t_land,
			t->drogue_validated_baro,
			t->main_validated_baro
	    );

	write_sd(&data_file, data_string);
}

void save_data_file(){
	flush_file(&data_file);
}

void write_datafile_message(char* msg){
	write_sd(&data_file, msg);
}

FRESULT write_headers(void)
{
    const char *header =
    	  "time,state,pressure,altitude,startAlt,temperature,"
          "velocity_world_x,velocity_world_y,velocity_world_z,"
    	  "baro_vz,lsm_accel_r,lsm_accel_p,lsm_accel_y,"
          "lsm_gyro_r,lsm_gyro_p,lsm_gyro_y,"
    	  "adxl_accel_r,adxl_accel_p,adxl_accel_y,"
    	  "predicted_apogee,airbrake_deployment,mag_r,mag_p,mag_y,"
    	  "q0,q1,q2,q3,accel_world_x,accel_world_y,accel_world_z,"
    	  "alt_fused,cam1_on,cam2_on,main_p_ematch_voltage,main_s_ematch_voltage,"
    	  "drogue_p_ematch_voltage,drogue_s_ematch_voltage,"
    	  "roll,pitch,yaw,t_burnout,t_apogee,t_drogue,t_main,t_land,"
    	  "drogue_validated_baro,main_validated_baro";

    return write_sd(&data_file, header);
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

uint32_t telemetry_log_period(FlightState_t state){
	switch(state) {
		case DISARMED: case LAUNCH_PAD: case LANDED:
			return 1000; // 1 Hz
		default:
			return 0; // As fast as possible
	}
}
