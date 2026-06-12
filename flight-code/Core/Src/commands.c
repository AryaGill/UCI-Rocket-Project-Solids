#include "commands.h"

#include "main.h"
#include "telemetry.h"
#include "sensors.h"
#include "cameras.h"
#include "buzzer.h"
#include "parachutes.h"
#include "airbrakes.h"
#include "telemetry.h"
#include <stdio.h>

char last_received_command[32];
uint32_t last_command_time = 0;

extern Bias_t bias;

//command handling from ground station
void handle_rf_command(char *cmd, FlightState_t *flight_state, Telemetry_t *telemetry) {
	// Check if receiving command from burst from gs
	uint32_t cur_time = HAL_GetTick();
	if (strcmp(last_received_command, cmd) == 0 && cur_time - last_command_time < 5000){
		return;
	}

	strcpy(last_received_command, cmd);
	last_command_time = cur_time;

	char line[64];
	snprintf(line, sizeof(line), "Received command: %s", cmd);
	snprintf(telemetry->cmd_echo, sizeof(telemetry->cmd_echo), cmd);
	write_datafile_message(line);
	//charge firing commands-> ONLY USED FOR TESTING
	if (strcmp(cmd, "Fire Main P") == 0){
	  	main_primary_on();
	  	HAL_Delay(CHARGE_DELAY);
	  	main_primary_off();
	} else if (strcmp(cmd, "Fire Main S") == 0){
		main_secondary_on();
		HAL_Delay(CHARGE_DELAY);
		main_secondary_off();
	} else if (strcmp(cmd, "Fire Drogue P") == 0){
		drogue_primary_on();
		HAL_Delay(CHARGE_DELAY);
		drogue_primary_off();
	} else if (strcmp(cmd, "Fire Drogue S") == 0){
		drogue_secondary_on();
		HAL_Delay(CHARGE_DELAY);
		drogue_secondary_off();

	//CAM control commands
	}  else if (strcmp(cmd, "CAM1ON") == 0){
		turn_camera_on(0);
	} else if (strcmp(cmd, "CAM2ON") == 0){
		turn_camera_on(1);
	} else if (strcmp(cmd, "CAM1OFF") == 0){
		turn_camera_off(0);
	} else if (strcmp(cmd, "CAM2OFF") == 0){
		turn_camera_off(1);
	} else if (strcmp(cmd, "ON") == 0){
			turn_camera_on(0);
			turn_camera_on(1);
	} else if (strcmp(cmd, "OFF") == 0){
			turn_camera_off(0);
			turn_camera_off(1);

	//Arms Rocket after it is on the rail-> USE BEFORE ALL LAUNCHES
	} else if (strcmp(cmd, "ARM") == 0){
		set_flight_state(LAUNCH_PAD, flight_state, telemetry);
		buzzer_set_frequency(4500);

	//plays servo airbrake sequence to verify extension
	} else if (strcmp(cmd, "SERVO SEQUENCE") == 0){
		perform_airbrakes_servo_sequence();

	//calibrate gyroscope bias once rocket is on the rail
	} else if (strcmp(cmd, "GYROCAL") == 0){
		Gyro_CalibrateBias(&bias, telemetry, 500);

	//clear SD card storage before launch
	} else if (strcmp(cmd, "RESET_SD") == 0){
		extern File_t data_file;
		sd_remove_file(&data_file);
		begin_data_file();
		snprintf(line, sizeof(line), "Received command: %s", cmd);
		snprintf(telemetry->cmd_echo, sizeof(telemetry->cmd_echo), cmd);
		write_datafile_message(line);
	}
}
