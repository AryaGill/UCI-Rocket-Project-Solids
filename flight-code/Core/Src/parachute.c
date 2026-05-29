/*
 * parachute.c
 *
 *  Created on: Mar 17, 2026
 *      Author: aryagill
 */


#include "parachute.h"
#include <math.h>

//verify parachute deployment was successful
void parachute_update_recovery_validation(FlightState_t *state, Telemetry_t *t){
	static uint16_t drogue_in_band_count = 0;
	static uint16_t main_in_band_count = 0;

	float vz = t->baro_vz;

	if (*state < DROGUE_PRIMARY_DEPLOYING){
		drogue_in_band_count = 0;
		main_in_band_count = 0;
		return;
	}

	// Ensure VZ is a valid value
	if (vz >= 0.0f || !isfinite(vz)){
		drogue_in_band_count = 0;
		main_in_band_count = 0;
		return;
	}

	// Drogue Validation
	if (!t->drogue_validated_baro){
		if (vz <= PARACHUTE_DROGUE_VEL_MAX && vz >= PARACHUTE_DROGUE_VEL_MIN){
			drogue_in_band_count++;
		} else {
			drogue_in_band_count = 0;
		}
		if (drogue_in_band_count >= PARACHUTE_DROGUE_SAMPLES_REQUIRED){
			t->drogue_validated_baro = 1;
		}
	}

	// Main Verification
	if (*state>=MAIN_PRIMARY_DEPLOYING && t->drogue_validated_baro == 1 && !t->main_validated_baro){
		if (vz <= PARACHUTE_MAIN_VEL_TARGET){
			main_in_band_count++;
		} else {
			main_in_band_count = 0;
		}
		if (main_in_band_count >= PARACHUTE_MAIN_SAMPLES_REQUIRED){
			t->main_validated_baro = 1;
		}
	}
}
