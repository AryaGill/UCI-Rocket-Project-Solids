#include "fsm.h"
#include "stm32h7xx_hal.h"
#include "sensors.h"
#include "main.h"
#include "telemetry.h"
#include "sd_card.h"
#include "telemetry.h"
#include "parachutes.h"
#include "airbrakes.h"
#include "launch_buffer.h"
#include <stdio.h>

// Some states need to know the time that the state started
uint32_t state_start_time = 0;

// Liftoff detection variables
uint32_t launch_accel_detected_time = -1;
unsigned int negative_accel_counter = 0;

// Motor Burn to Gliding Ascent detection variables
uint32_t num_increasing_accel = 0;
float prev_accel = 0;

uint32_t prev_led_toggle_time = 0;

extern launch_buffer_t launch_buffer;
extern File_t data_file;

void set_flight_state(FlightState_t new_state, FlightState_t *flight_state, Telemetry_t *telemetry) {
	*flight_state = new_state;

	// Write
	char line[50];
	char state_str[30];
	state_to_string_name(new_state, state_str);
	snprintf(line, sizeof(line), "Entering state: %s", state_str);
	write_datafile_message(line);
	write_sd_state(STATE_FILE, *flight_state, telemetry->startAlt);
}

uint8_t sensors_indicate_flight(Telemetry_t *telemetry){
	float alt_i = telemetry->altitude;

	for (int i = 0; i < 250; ++i){
		read_sensors(telemetry);
		if (abs(telemetry->altitude - alt_i) > POWER_RESET_MIN_ALT_CHANGE){
			return 1;
		}
		HAL_Delay(10);
	}

	return 0;
}

void init_flight_state(FlightState_t *flight_state, Telemetry_t *telemetry) {
	telemetry->alt_fused = 0;
	
	//possible start altitude after reset fix
	for (int i = 0; i < 10; ++i) {
		LPS22HH_Read(telemetry);
		if (telemetry->temperature != -999) {
			telemetry->startAlt = telemetry->altitude;
		}
		HAL_Delay(10);
	}

	// Detect power reset
	if (sd_file_exists(STATE_FILE)){
		FlightState_t sd_state;
		float sd_start_alt;
		read_sd_state(STATE_FILE, &sd_state, &sd_start_alt);

		uint32_t reset_flags = RCC->RSR;
		// If Power-on reset (power removed and restored) and sensors indicate in flight and alt > threshold
		if ((reset_flags & RCC_RSR_PORRSTF) && (telemetry->altitude - sd_start_alt > MIN_RESET_ALT) && sensors_indicate_flight(telemetry) == 1){
			// Print message in data file
			write_datafile_message("POWER RESET DETECTED");

			// Go to correct state and start altitude
			set_flight_state(sd_state, flight_state, telemetry);
			telemetry->startAlt = sd_start_alt;

			return;
		}
	}

	// Remove Reset Flags
	RCC->RSR |= RCC_RSR_RMVF;

	set_flight_state(DISARMED, flight_state, telemetry);
}

void update_flight_state(FlightState_t *flight_state, Telemetry_t *telemetry) {
	// Determine Next State
	switch(*flight_state) {
		case DISARMED:
			// Toggle LED
			uint32_t cur_time = HAL_GetTick();
			if (cur_time - prev_led_toggle_time > 100){
				HAL_GPIO_TogglePin(LED_GPIO_Port, LED_Pin);
				prev_led_toggle_time = cur_time;
			}
			break;
    	case LAUNCH_PAD:
    		// Detect if launched

    		if (launch_accel_detected_time == -1){
				// Acceleration not detected yet
				//changed to < for vacuum testing
				if (telemetry->bmx_accel_r > LAUNCH_ACCEL_THRESHOLD){
					// Positive acceleration detected. Begin period of waiting to get off rail.
					launch_accel_detected_time = HAL_GetTick();
					negative_accel_counter = 0;
				}
    		}
			else if (HAL_GetTick() - launch_accel_detected_time > RAIL_DELAY_TIME){
				// In evaluation period. Monitor for any negative acceleration value.
				// If detected, reset the system and begin again.
				if (HAL_GetTick() - launch_accel_detected_time > RAIL_DELAY_TIME + LAUNCH_EVAL_PERIOD_TIME){
					// Enough time passed without negative acceleration. Launch detected
					set_flight_state(MOTOR_BURN, flight_state, telemetry);
					write_sd(&data_file, "launch buffer dump");
					launch_buffer_flush(&launch_buffer);
					state_start_time = HAL_GetTick();
				}
				else if (telemetry->bmx_accel_r < 9.81){ // 0
					// Negative acceleration detected. Reset system.
					if (++negative_accel_counter >= 5){ // 20
						launch_accel_detected_time = -1;
						write_datafile_message("LAUNCH DETECTION FAILED");
					}
				}
				else {
					negative_accel_counter = 0;
				}
			}

    		break;

    	case MOTOR_BURN:
    		if (telemetry->accel_world_z > prev_accel){
    			++num_increasing_accel;
    		}
    		else{
    			num_increasing_accel = 0;
    		}
    		prev_accel = telemetry->accel_world_z;

    		if (HAL_GetTick() - state_start_time > MOTOR_BURN_TIME) {
    			set_flight_state(GLIDING_ASCENT, flight_state, telemetry);
    		}
    		else if (num_increasing_accel > 10 && telemetry->accel_world_z < 0) {
    			write_datafile_message("MOTOR BURNOUT DETECTED");
    			set_flight_state(GLIDING_ASCENT, flight_state, telemetry);
    		}

    		break;

    	case GLIDING_ASCENT:
    		if (telemetry->baro_vz < APOGEE_VELO_THRESHOLD && telemetry->altitude - telemetry->startAlt > DROGUE_DEPLOY_MIN_ALT) {
    			set_flight_state(DROGUE_PRIMARY_DEPLOYING, flight_state, telemetry);
    			drogue_primary_on();
    			//drogue primary starts firing and the time this starts is stored
    			state_start_time = HAL_GetTick();
    		}

    		break;

    	case DROGUE_PRIMARY_DEPLOYING:
    		if (HAL_GetTick() - state_start_time >= CHARGE_DELAY) {
    			drogue_primary_off();

    			//time that primary finishes is stored and bool is set to true so this state does not run again
    			state_start_time = HAL_GetTick();

    			set_flight_state(DROGUE_PRIMARY_DEPLOYED, flight_state, telemetry);
    		}

    		break;

    	case DROGUE_PRIMARY_DEPLOYED:
    		if (HAL_GetTick() - state_start_time >= BACKUP_DELAY) {
    			drogue_secondary_on();

    			//time when secondary finishes is stored and bools set so this state does not run again
    			state_start_time = HAL_GetTick();

    			set_flight_state(DROGUE_SECONDARY_DEPLOYING, flight_state, telemetry);
    		}

    		break;

    	case DROGUE_SECONDARY_DEPLOYING:
    		if (HAL_GetTick() - state_start_time >= CHARGE_DELAY){
    			drogue_secondary_off();

    			set_flight_state(DROGUE_SECONDARY_DEPLOYED, flight_state, telemetry);
    		}

    		break;

    	case DROGUE_SECONDARY_DEPLOYED:
    		// Wait for main deployment
    		if (telemetry->altitude - telemetry->startAlt < MAIN_DEPLOY_MAX_ALT && telemetry->altitude - telemetry->startAlt > MAIN_DEPLOY_MIN_ALT){
    			main_primary_on();
    			state_start_time = HAL_GetTick();
    			set_flight_state(MAIN_PRIMARY_DEPLOYING, flight_state, telemetry);
    		}

    		break;

    	case MAIN_PRIMARY_DEPLOYING:
    		if(HAL_GetTick() - state_start_time >= CHARGE_DELAY){
    			main_primary_off();
    			state_start_time = HAL_GetTick();
    			set_flight_state(MAIN_PRIMARY_DEPLOYED, flight_state, telemetry);
    		}

    		break;

    	case MAIN_PRIMARY_DEPLOYED:
    		if(HAL_GetTick() - state_start_time >= BACKUP_DELAY){
    			main_secondary_on();
    			state_start_time = HAL_GetTick();
    			set_flight_state(MAIN_SECONDARY_DEPLOYING, flight_state, telemetry);
    		}

    		break;

    	case MAIN_SECONDARY_DEPLOYING:
    		if(HAL_GetTick() - state_start_time >= CHARGE_DELAY){
    			main_secondary_off();
    			set_flight_state(MAIN_SECONDARY_DEPLOYED, flight_state, telemetry);
    		}
    		break;

    	case MAIN_SECONDARY_DEPLOYED:
    		if (telemetry->baro_vz > LANDED_VELO_THRESHOLD){ // Change condition // TODO Needs to be switched according to chat?
    			set_flight_state(LANDED, flight_state, telemetry);
    			sd_delete_file(STATE_FILE);
    		}
    		break;

    	case LANDED:

    		break;
	}
}
