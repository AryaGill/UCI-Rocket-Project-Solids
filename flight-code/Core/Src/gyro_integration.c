#include "gyro_integration.h"

Quaternions_t quats[2];
uint8_t active_quats = 0;
uint32_t last_quat_change_time;

void Gyro_Integrate(Telemetry_t* telemetry, Quaternions_t* q)
{
    static uint64_t prev_time = 0;

    uint64_t now = micros();
    float dt = (now - prev_time) * 1e-6f;
    prev_time = now;

    if (dt <= 0.0f || dt > 0.02f)
        return;

    float q0 = quats->q0;
    float q1 = quats->q1;
    float q2 = quats->q2;
    float q3 = quats->q3;

    // Gyro rates (rad/s)
    float gx = telemetry->lsm_gyro_p;
    float gy = telemetry->lsm_gyro_y;
    float gz = telemetry->lsm_gyro_r;

    // Quaternion derivative
    float qDot0 = 0.5f * (-q1*gx - q2*gy - q3*gz);
    float qDot1 = 0.5f * ( q0*gx + q2*gz - q3*gy);
    float qDot2 = 0.5f * ( q0*gy - q1*gz + q3*gx);
    float qDot3 = 0.5f * ( q0*gz + q1*gy - q2*gx);

    // Integrate
    q0 += qDot0 * dt;
    q1 += qDot1 * dt;
    q2 += qDot2 * dt;
    q3 += qDot3 * dt;

    // Normalize quaternion
    float norm = sqrtf(q0*q0 + q1*q1 + q2*q2 + q3*q3);
    float inv = 1.0f / norm;

    q->q0 = q0 * inv;
    q->q1 = q1 * inv;
    q->q2 = q2 * inv;
    q->q3 = q3 * inv;
}

void Quaternion_From_Accel(Telemetry_t* telemetry, Quaternions_t* q)
{
    float ax = telemetry->lsm_accel_p;
    float ay = telemetry->lsm_accel_y;
    float az = telemetry->lsm_accel_r;

    // Normalize acceleration
    float norm = sqrtf(ax*ax + ay*ay + az*az);
    if (norm == 0.0f) return;

    ax /= norm;
    ay /= norm;
    az /= norm;

    // Gravity in world frame (Z-up)
    float gx = 0.0f;
    float gy = 0.0f;
    float gz = 1.0f;

    // Dot product (cos of angle between vectors)
    float dot = gx*ax + gy*ay + gz*az;

    float q0, q1, q2, q3;

    if (dot < -0.999999f)
    {
        // 180° rotation: choose arbitrary axis perpendicular to gravity
        q0 = 0.0f;
        q1 = 1.0f;
        q2 = 0.0f;
        q3 = 0.0f;
    }
    else
    {
        float s = sqrtf((1.0f + dot) * 2.0f);
        float invs = 1.0f / s;

        // Cross product g × a
        float cx = gy*az - gz*ay;
        float cy = gz*ax - gx*az;
        float cz = gx*ay - gy*ax;

        q0 = s * 0.5f;
        q1 = cx * invs;
        q2 = cy * invs;
        q3 = cz * invs;
    }

    q->q0 = q0;
    q->q1 = q1;
    q->q2 = q2;
    q->q3 = q3;
}

void integrate_gyro(FlightState_t flight_state, Telemetry_t* telemetry){
	switch (flight_state){
	case LAUNCH_PAD: case DISARMED:
		// Switch active set of quats after 5 seconds
		uint32_t cur_time = HAL_GetTick();
		if (cur_time - last_quat_change_time > 5000){
			last_quat_change_time = cur_time;
			active_quats = (active_quats == 0) ? 1 : 0;
		}

		// Integrate gyro for active quats
		Gyro_Integrate(telemetry, &quats[active_quats]);

		// Get quats from accel on inactive quats
		uint8_t inactive_quats = (active_quats == 0) ? 1 : 0;
		Quaternion_From_Accel(telemetry, &quats[inactive_quats]);
	default:
		// Integrate gyro for just active set of quats
		Gyro_Integrate(telemetry, &quats[active_quats]);
	}

	// Set quats to active set of quats
	telemetry->q0 = quats[active_quats].q0;
	telemetry->q1 = quats[active_quats].q1;
	telemetry->q2 = quats[active_quats].q2;
	telemetry->q3 = quats[active_quats].q3;
}
