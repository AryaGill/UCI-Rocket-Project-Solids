#include "complementary_filter.h"
#include "airbrakes.h"
#include <math.h>
#include <stdint.h>
#include <stdbool.h>

static uint64_t prev_time_cf_us = 0;

static float prev_baro_alt = 0.0f;
static float baro_velocity_filt = 0.0f;

static bool cf_initialized = false;

//Initializes Complementary filter values
void complementary_filter_init(Telemetry_t *telemetry)
{
    float baro_alt = telemetry->altitude - telemetry->startAlt;

    prev_time_cf_us = micros();
    prev_baro_alt = baro_alt;
    baro_velocity_filt = 0.0f;

    telemetry->velocity_world_x = 0.0f;
    telemetry->velocity_world_y = 0.0f;
    telemetry->velocity_world_z = 0.0f;

    telemetry->baro_vz = 0.0f;
    telemetry->alt_fused = baro_alt;

    telemetry->time_until_trust_baro = 0.0f;
    telemetry->prev_deployment = telemetry->airbrake_deployment;

    cf_initialized = true;
}

/**
 * Updates fused altitude and vertical velocity estimates.
 *
 * telemetry: flight telemetry struct to update
 * flight_state: current flight state
 *
 * Combines IMU acceleration with barometer altitude/velocity.
 */
void complementary_filter(Telemetry_t *telemetry, FlightState_t *flight_state)
{
    if (!cf_initialized) {
        complementary_filter_init(telemetry);
        return;
    }

    //calculate change in time (time step)
    uint64_t cur_time_us = micros();
    float dt = (float)(cur_time_us - prev_time_cf_us) * 1e-6f;
    prev_time_cf_us = cur_time_us;

    //skips computation if time step is invalid
    if (!isfinite(dt) || dt <= 0.0f || dt > 0.1f) {
        return;
    }

    // Determine weight of barometer measurements based on time since last airbrakes change
    if (telemetry->airbrake_deployment == telemetry->prev_deployment){
    	telemetry->time_until_trust_baro = fmaxf(0.0f, telemetry->time_until_trust_baro - dt);
    }else {
    	telemetry->time_until_trust_baro += BARO_TRUST_TIME * fabs(telemetry->prev_deployment - telemetry->airbrake_deployment) / (NUM_DEPLOYMENT_LEVELS - 1);
        telemetry->time_until_trust_baro = fminf(BARO_TRUST_TIME, telemetry->time_until_trust_baro);
    }
    telemetry->prev_deployment = telemetry->airbrake_deployment;

    float x = fmaxf(0.0f, fminf(1.0f, 1.0f - (telemetry->time_until_trust_baro / BARO_TRUST_TIME)));
    float w_baro = x * x * (3 - 2 * x);

    // Calculate barometric velocity
    float baro_alt = telemetry->altitude - telemetry->startAlt;

    float velocity_baro_raw = (baro_alt - prev_baro_alt) / dt;
    prev_baro_alt = baro_alt;

    float baro_vel_alpha = TAU_BARO_VEL / (TAU_BARO_VEL + dt);
    telemetry->baro_vz =
        baro_vel_alpha * telemetry->baro_vz +
        (1.0f - baro_vel_alpha) * velocity_baro_raw;

    // Calculate velocity_world_z
    float alpha_velocity = TAU_VELOCITY / (TAU_VELOCITY + dt);
	float velocity_imu = telemetry->velocity_world_z + telemetry->accel_world_z * dt;

	float velocity_fused =
	    	alpha_velocity * velocity_imu +
	        (1.0f - alpha_velocity) * telemetry->baro_vz;

	telemetry->velocity_world_z = (1-w_baro) * velocity_imu + w_baro * velocity_fused;

    // Calculate fused altitude
    float altitude_pred =
        telemetry->alt_fused + telemetry->velocity_world_z * dt;

    float alpha_altitude = TAU_ALTITUDE / (TAU_ALTITUDE + dt);
    float altitude_fused =
    	alpha_altitude * altitude_pred +
        (1.0f - alpha_altitude) * baro_alt;

    telemetry->alt_fused = (1-w_baro) * altitude_pred + w_baro * altitude_fused;
    if (!isfinite(telemetry->velocity_world_z)) telemetry->velocity_world_z = 0.0f;
    if (!isfinite(telemetry->alt_fused))        telemetry->alt_fused = baro_alt;
}
