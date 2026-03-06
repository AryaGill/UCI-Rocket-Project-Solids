#include "commands.h"

#include "main.h"
#include "telemetry.h"
#include "sensors.h"
#include "cameras.h"
#include "buzzer.h"
#include "parachutes.h"
#include "airbrakes.h"
#include <stdio.h>

char last_received_command[32];
uint32_t last_command_time = 0;

void handle_rf_command(char *cmd, FlightState_t *flight_state, Telemetry_t *telemetry) {
	// Check if receiving command from burst from gs
	uint32_t cur_time = HAL_GetTick();
	if (strcmp(last_received_command, cmd) == 0 && cur_time - last_command_time < 2300){
		return;
	}

	strcpy(last_received_command, cmd);
	last_command_time = cur_time;

	char line[64];
	snprintf(line, sizeof(line), "Received command: %s", cmd);
	write_sd(FLIGHT_DATA_FILE, line);
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
	}  else if (strcmp(cmd, "CAM1ON") == 0){
		turn_camera_on(0);
	} else if (strcmp(cmd, "CAM2ON") == 0){
		turn_camera_on(1);
	} else if (strcmp(cmd, "CAM1OFF") == 0){
		turn_camera_off(0);
	} else if (strcmp(cmd, "CAM2OFF") == 0){
		turn_camera_off(1);
	} else if (strcmp(cmd, "ARM") == 0){
		set_flight_state(LAUNCH_PAD, flight_state, telemetry);
		buzzer_set_frequency(4500);
	}
	else if (strcmp(cmd, "SERVO SEQUENCE") == 0){
		perform_airbrakes_servo_sequence();
	}
}
