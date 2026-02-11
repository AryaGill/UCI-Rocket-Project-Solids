#include "commands.h"

#include "main.h"
#include "telemetry.h"
#include "sensors.h"
#include <stdio.h>

void handle_rf_command(char *cmd, FlightState_t *flight_state, Telemetry_t *telemetry) {
	char line[30];
	snprintf(line, sizeof(line), "Received command: %s", cmd);
	write_sd(FLIGHT_DATA_FILE, line);
	if(strcmp(cmd, "ON") == 0){
//		Serial.println("Camera On Recieved");
//	    HWSERIAL.println("TEENSY Camera on");
//		digitalWrite(camera1,HIGH);
//		digitalWrite(camera2,HIGH);
	}else if (strcmp(cmd, "OFF") == 0){
//		Serial.println("Camera Off Recieved");
//		HWSERIAL.println("TEENSY Camera off");
//		digitalWrite(camera1, LOW);
//		digitalWrite(camera2, LOW);
	}else if (strcmp(cmd, "Fire Main P") == 0){
	  	HAL_GPIO_WritePin(Main_Parachute_1_GPIO_Port, Main_Parachute_1_Pin, GPIO_PIN_SET);
	  	HAL_Delay(CHARGE_DELAY);
	  	HAL_GPIO_WritePin(Main_Parachute_1_GPIO_Port, Main_Parachute_1_Pin, GPIO_PIN_RESET);
	} else if (strcmp(cmd, "Fire Main S") == 0){
		HAL_GPIO_WritePin(Main_Parachute_2_GPIO_Port, Main_Parachute_2_Pin, GPIO_PIN_SET);
		HAL_Delay(CHARGE_DELAY);
		HAL_GPIO_WritePin(Main_Parachute_2_GPIO_Port, Main_Parachute_2_Pin, GPIO_PIN_RESET);
	} else if (strcmp(cmd, "Fire Drogue P") == 0){
		HAL_GPIO_WritePin(Drogue_Parachute_1_GPIO_Port, Drogue_Parachute_1_Pin, GPIO_PIN_SET);
		HAL_Delay(CHARGE_DELAY);
		HAL_GPIO_WritePin(Drogue_Parachute_1_GPIO_Port, Drogue_Parachute_1_Pin, GPIO_PIN_RESET);
	} else if (strcmp(cmd, "Fire Drogue S") == 0){
		HAL_GPIO_WritePin(Drogue_Parachute_2_GPIO_Port, Drogue_Parachute_2_Pin, GPIO_PIN_SET);
		HAL_Delay(CHARGE_DELAY);
		HAL_GPIO_WritePin(Drogue_Parachute_2_GPIO_Port, Drogue_Parachute_2_Pin, GPIO_PIN_RESET);
	}  else if (strcmp(cmd, "CAM1ON") == 0){
		turn_camera_on(0);
	} else if (strcmp(cmd, "CAM2ON") == 0){
		turn_camera_on(1);
	} else if (strcmp(cmd, "CAM1OFF") == 0){
		turn_camera_off(0);
	} else if (strcmp(cmd, "CAM2OFF") == 0){
		turn_camera_on(1);
	} else if (strcmp(cmd, "ARM") == 0){
		set_flight_state(LAUNCH_PAD, flight_state, telemetry);
	}
}
