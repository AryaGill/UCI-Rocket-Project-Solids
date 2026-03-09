#include "complementary_filter.h"

uint32_t prev_time_cf = 0;

float prev_baro_alt = 0.0f;  // previous barometer altitude for velocity calculation

void update_horizontal_velocity(Telemetry_t* telemetry, float dt) {
    // Earth-frame horizontal acceleration (tilt-compensated)
    float ax = telemetry->accel_world_x;
    float ay = telemetry->accel_world_y;

    // Time-based decay constant
    // alpha_per_second = fraction of velocity retained per second
    const float alpha_per_second = 0.95f; // adjust: 0.98 → slow decay, 0.95 → faster decay

    // Convert to per-update decay
    float alpha_dt = powf(alpha_per_second, dt);

    // Update velocities
    telemetry->velocity_world_x = alpha_dt * telemetry->velocity_world_x + (1.0f - alpha_dt) * ax * dt;
    telemetry->velocity_world_y = alpha_dt * telemetry->velocity_world_y + (1.0f - alpha_dt) * ay * dt;
}

void complementary_filter(Telemetry_t* telemetry) {
	// Get time since last call
    uint32_t cur_time = micros();
    float dt = (cur_time - prev_time_cf) * 1e-6f;
    if (dt <= 0.0f) return;  // safety check
    if (dt > 0.05f) dt = 0.05f; // Clamp large dt
    prev_time_cf = cur_time;

    // STAGE 1: Velocity Fusion
    // Use gravity-compensated, tilt-corrected vertical acceleration
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
    float alt_from_velocity = telemetry->alt_fused + telemetry->velocity_world_z * dt;

    // Fuse altitudes: 95% integrated velocity, 5% raw barometer
    telemetry->alt_fused = ALPHA_ALTITUDE * alt_from_velocity
			  + (1.0f - ALPHA_ALTITUDE) * baro_alt;

   update_horizontal_velocity(telemetry, dt);
}
