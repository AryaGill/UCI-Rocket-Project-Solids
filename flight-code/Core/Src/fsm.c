#include "fsm.h"
#include "stm32h7xx_hal.h"
#include "telemetry.h"
#include "sensors.h"
#include "main.h"

float alt_dif_buffer[ALT_DIF_BUF_SIZE];
int alt_dif_buffer_idx = 0;
int prev_alt_time = 0;
float prev_alt;

// Time variables
unsigned long launch_start_time = 0;
unsigned long drogue_primary_start_time = 0;
unsigned long drogue_primary_end_time = 0;
unsigned long drogue_secondary_start_time = 0;
unsigned long main_primary_start_time = 0;
unsigned long main_primary_end_time = 0;
unsigned long main_secondary_start_time = 0;

float get_avg_alt_dif() {
	float sum = 0;
	float largest = alt_dif_buffer[0];
	float smallest = alt_dif_buffer[0];
	for (int i = 0; i < ALT_DIF_BUF_SIZE; ++i){
		sum += alt_dif_buffer[i];
		largest = max(largest, alt_dif_buffer[i]);
    	smallest = min (smallest, alt_dif_buffer[i]);
	}
  	return (sum - largest - smallest) / (ALT_DIF_BUF_SIZE - 2);
}

void update_alt_dif_buf(float new_alt_dif) {

	float cur_time = HAL_GetTick();
	if (cur_time <= prev_alt_time){
		return;
	}
	alt_dif_buffer[alt_dif_buffer_idx] = new_alt_dif / (cur_time - prev_alt_time) * 1000;
	alt_dif_buffer_idx = (alt_dif_buffer_idx + 1) % ALT_DIF_BUF_SIZE;
	prev_alt_time = cur_time;
}

void set_flight_state(FlightState_t new_state, FlightState_t *flight_state) {
	*flight_state = new_state;

	// Write flight state to state file
//	stateFile = SD.open("rocket_state.csv", FILE_WRITE);
//	if (stateFile) {
//		stateFile.println(String(startAlt, 8));
//		stateFile.flush();
//      	stateFile.println((int)flight_state);
//      	stateFile.flush();
//	} else {
//		Serial.println("Failed to create state file");
//	}
}

void initialize_flight_state(FlightState_t *flight_state, Telemetry_t *telemetry) {
	//possible start altitude after reset fix
	for (int i = 0; i < ALT_DIF_BUF_SIZE + 1; ++i) {
		read_bmp(telemetry);
		if (telemetry->temperature != -999) {
			telemetry->startAlt = telemetry->altitude;
			update_alt_dif_buf(telemetry->startAlt - prev_alt);
			prev_alt = telemetry->startAlt;
		}
		else {
//			Serial.println("Failed to get BPM390 data");
		}
	}

	// Open rocket state file
//	stateFile = SD.open("rocket_state.csv", FILE_READ);
//
//	if (stateFile) {
//		if (stateFile.size() > 0) {
//			float read_alt = startAlt;
//			// File exists and has data
//			startAlt = stateFile.readStringUntil('\n').trim().toFloat();
//			Serial.print("startAlt loaded from file: ");
//			Serial.println(startAlt);
//
//			if (read_alt - startAlt > 183){
//				flight_state = static_cast<FlightState>(stateFile.readStringUntil('\n').trim().toInt());
//				Serial.print("Flight state loaded from file: ");
//				Serial.println((int)flight_state);
//			}
//			else{
//				set_flight_state(LAUNCH_PAD);
//			}
//
//		} else {
//			// File exists but empty, close previous read mode
//			Serial.println("Writing starting altitude.");
//			stateFile.close();
//			set_flight_state(LAUNCH_PAD);
//		}
//	} else {
//		 //File doesn't exist
//		Serial.println("Creating and writing starting altitude.");
		set_flight_state(LAUNCH_PAD, flight_state);
//	}
}

void update_flight_state(FlightState_t *flight_state, Telemetry_t *telemetry) {
	update_alt_dif_buf(telemetry->altitude - prev_alt);

	// Determine Next State
	switch(*flight_state) {
    	case LAUNCH_PAD:
    		// Detect if launched
    		if (get_avg_alt_dif() > LAUNCH_THRESHOLD) {
//    			dataFile.println("LAUNCHED");
    			launch_start_time = HAL_GetTick();
    			set_flight_state(MOTOR_BURN, flight_state);
    		}

    		break;

    	case MOTOR_BURN:
    		// Add logic: wait for certain delay or for acceleration to change
    		if (HAL_GetTick() - launch_start_time > 1500){
    			set_flight_state(GLIDING_ASCENT, flight_state);
    		}

    		break;

    	case GLIDING_ASCENT:
    		if (get_avg_alt_dif() < APOGEE_THRESHOLD) {
    			set_flight_state(DROGUE_PRIMARY_DEPLOYING, flight_state);
    			HAL_GPIO_WritePin(Drogue_Parachute_1_GPIO_Port, Drogue_Parachute_1_Pin, GPIO_PIN_SET);
//    			dataFile.println("Primary Drogue Deployed");
    			//drogue primary starts firing and the time this starts is stored
    			drogue_primary_start_time = HAL_GetTick();
    		}

    		break;

    	case DROGUE_PRIMARY_DEPLOYING:
    		if (HAL_GetTick() - drogue_primary_start_time >= CHARGE_DELAY) {
    			HAL_GPIO_WritePin(Drogue_Parachute_1_GPIO_Port, Drogue_Parachute_1_Pin, GPIO_PIN_RESET);

    			//time that primary finishes is stored and bool is set to true so this state does not run again
    			drogue_primary_end_time = HAL_GetTick();

    			set_flight_state(DROGUE_PRIMARY_DEPLOYED, flight_state);
    		}

    		break;

    	case DROGUE_PRIMARY_DEPLOYED:
    		if (HAL_GetTick() - drogue_primary_end_time >= backup_delay) {
    			HAL_GPIO_WritePin(Drogue_Parachute_2_GPIO_Port, Drogue_Parachute_2_Pin, GPIO_PIN_SET);
//    			dataFile.println("Secondary Drogue Deployed");

    			//time when secondary finishes is stored and bools set so this state does not run again
    			drogue_secondary_start_time = HAL_GetTick();

    			set_flight_state(DROGUE_SECONDARY_DEPLOYING, flight_state);
    		}

    		break;

    	case DROGUE_SECONDARY_DEPLOYING:
    		if (HAL_GetTick() - drogue_secondary_start_time >= CHARGE_DELAY){
    			HAL_GPIO_WritePin(Drogue_Parachute_2_GPIO_Port, Drogue_Parachute_2_Pin, GPIO_PIN_RESET);

    			set_flight_state(DROGUE_SECONDARY_DEPLOYED, flight_state);
    		}

    	case DROGUE_SECONDARY_DEPLOYED:
    		// Wait for main deployment
    		if (telemetry->altitude - telemetry->startAlt < 229 && telemetry->altitude - telemetry->startAlt > 77){
    			HAL_GPIO_WritePin(Main_Parachute_1_GPIO_Port, Main_Parachute_1_Pin, GPIO_PIN_SET);
//    			dataFile.println("Primary Main Deployed");
    			main_primary_start_time = HAL_GetTick();
    			set_flight_state(MAIN_PRIMARY_DEPLOYING, flight_state);
    		}

    		break;

    	case MAIN_PRIMARY_DEPLOYING:
    		if(millis() - main_primary_start_time >= CHARGE_DELAY){
    			HAL_GPIO_WritePin(Main_Parachute_1_GPIO_Port, Main_Parachute_1_Pin, GPIO_PIN_RESET);
    			main_primary_end_time = HAL_GetTick();
    			set_flight_state(MAIN_PRIMARY_DEPLOYED, flight_state);
    		}

    		break;

    	case MAIN_PRIMARY_DEPLOYED:
    		if(HAL_GetTick() - main_primary_end_time >= BACKUP_DELAY){
    			HAL_GPIO_WritePin(Main_Parachute_2_GPIO_Port, Main_Parachute_2_Pin, GPIO_PIN_SET);
//    			dataFile.println("Secondary Main Deployed");
    			main_secondary_start_time = HAL_GetTick();
    			set_flight_state(MAIN_SECONDARY_DEPLOYING, flight_state);
    		}

    		break;

    	case MAIN_SECONDARY_DEPLOYING:
    		if(HAL_GetTick() - main_secondary_start_time >= CHARGE_DELAY){
    			HAL_GPIO_WritePin(Main_Parachute_2_GPIO_Port, Main_Parachute_2_Pin, GPIO_PIN_RESET);
    			set_flight_state(MAIN_SECONDARY_DEPLOYED, flight_state);
    		}
    		break;

    	case MAIN_SECONDARY_DEPLOYED:
    		if (get_avg_alt_dif() > LANDED_THRESHOLD){ // Change condition
    			set_flight_state(LANDED, flight_state);
    		}
    		break;

    	case LANDED:

    		break;
	}

	prev_alt = telemetry->altitude;
}

String state_to_string(FlightState state) {
	switch(state) {
    	case LAUNCH_PAD:
    		return "LAUNCH_PAD";
    		break;
    	case MOTOR_BURN:
    		return "MOTOR_BURN";
    		break;
    	case GLIDING_ASCENT:
    		return "GLIDING_ASCENT";
    		break;
    	case DROGUE_PRIMARY_DEPLOYING:
    		return "DROGUE_PRIMARY_DEPLOYING";
    		break;
    	case DROGUE_PRIMARY_DEPLOYED:
    		return "DROGUE_PRIMARY_DEPLOYED";
    		break;
    	case DROGUE_SECONDARY_DEPLOYING:
    		return "DROGUE_SECONDARY_DEPLOYING";
    		break;
    	case DROGUE_SECONDARY_DEPLOYED:
    		return "DROGUE_SECONDARY_DEPLOYED";
    		break;
    	case MAIN_PRIMARY_DEPLOYING:
    		return "MAIN_PRIMARY_DEPLOYING";
    		break;
    	case MAIN_PRIMARY_DEPLOYED:
    		return "MAIN_PRIMARY_DEPLOYED";
    		break;
    	case MAIN_SECONDARY_DEPLOYING:
    		return "MAIN_SECONDARY_DEPLOYING";
    		break;
    	case MAIN_SECONDARY_DEPLOYED:
    		return "MAIN_SECONDARY_DEPLOYED";
    		break;
    	case LANDED:
    		return "LANDED";
    		break;
    	default:
    		return "UNKNOWN_STATE";
    		break;
	}
}
