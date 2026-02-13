#include "complementary_filter.h"

uint32_t prev_time_cf = 0;

float alt_fused = 0.0f;      // fused altitude (m)
float prev_baro_alt = 0.0f;  // previous barometer altitude for velocity calculation

void complementary_filter(Telemetry_t* telemetry) {
	// Get time since last call
    uint32_t cur_time = micros();
    float dt = (cur_time - prev_time_cf) * 1e-6f;
    if (dt <= 0.0f) return;  // safety check
    if (dt > 0.05f) dt = 0.05f; // Clamp large dt
    prev_time_cf = cur_time;

    // STAGE 1: Velocity Fusion
    // Use gravity-compensated, tilt-corrected vertical acceleration
    telemetry->velocity_world_x += telemetry->accel_world_x * dt;
    telemetry->velocity_world_y += telemetry->accel_world_y * dt;
    telemetry->velocity_world_z += telemetry->accel_world_z * dt;
    // float velocity_imu = velocity_fused + (-Accel_z - 9.81f) * dt;

    // Calculate barometric velocity
    float baro_alt = telemetry->altitude - telemetry->startAlt;
    float velocity_baro = (baro_alt - prev_baro_alt) / dt;
    prev_baro_alt = baro_alt;

    // Fuse velocities: 99% IMU, 1% barometer
    telemetry->velocity_world_z = ALPHA_VELOCITY * telemetry->velocity_world_z
				 + (1.0f - ALPHA_VELOCITY) * velocity_baro;

    // STAGE 2: Altitude Fusion
    // Integrate fused velocity to get altitude prediction
    float alt_from_velocity = alt_fused + telemetry->velocity_world_z * dt;

    // Fuse altitudes: 95% integrated velocity, 5% raw barometer
    telemetry->alt_fused = ALPHA_ALTITUDE * alt_from_velocity
			  + (1.0f - ALPHA_ALTITUDE) * baro_alt;

    // Damp x and y velocity to reduce drift
    telemetry->velocity_world_x *= 0.999f;
    telemetry->velocity_world_y *= 0.999f;
}
