#include "launch_buffer.h"
#include "sd_card.h"
#include "telemetry.h"
#include <stdio.h>
#include <string.h>

extern File_t data_file;

void launch_buffer_init(launch_buffer_t *lb)
{
    lb->head = 0;
    lb->full = 0;
    lb->flushed = 0;
}

void launch_buffer_add(launch_buffer_t *lb, FlightState_t flight_state, Telemetry_t *t)
{
	if (lb->flushed) return;
    char *data_string = lb->buffer[lb->head];
    char state_str[3];

    state_to_string_num(flight_state, state_str);

    snprintf(data_string, LAUNCH_BUFFER_LINE_SIZE,
        "%lu,%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%u,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%u,%u,%lu,%lu,%lu,%lu,%.3f,%.3f,%.3f,%lu,%lu,%lu,%lu,%lu,%i,%i",
        t->time,
        state_str,
        t->pressure,
        t->altitude,
        t->startAlt,
        t->temperature,
        t->velocity_world_x,
        t->velocity_world_y,
        t->velocity_world_z,
        t->baro_vz,
        t->lsm_accel_r,
        t->lsm_accel_p,
        t->lsm_accel_y,
        t->lsm_gyro_r,
        t->lsm_gyro_p,
        t->lsm_gyro_y,
        t->adxl_accel_r,
        t->adxl_accel_p,
        t->adxl_accel_y,
        t->predicted_apogee,
        t->airbrake_deployment,
        t->mag_r,
        t->mag_p,
        t->mag_y,
        t->q0,
        t->q1,
        t->q2,
        t->q3,
        t->accel_world_x,
        t->accel_world_y,
        t->accel_world_z,
        t->alt_fused,
        t->cam1_on,
        t->cam2_on,
        t->main_p_ematch_voltage,
        t->main_s_ematch_voltage,
        t->drogue_p_ematch_voltage,
        t->drogue_s_ematch_voltage,
        t->roll,
        t->pitch,
        t->yaw,
        t->t_burnout,
        t->t_apogee,
        t->t_drogue,
        t->t_main,
        t->t_land,
        t->drogue_validated_baro,
        t->main_validated_baro
    );

    lb->head++;

    if (lb->head >= LAUNCH_BUFFER_SIZE) {
        lb->head = 0;
        lb->full = 1;
    }
}

void launch_buffer_flush(launch_buffer_t *lb)
{
    if (lb->flushed) return;

    uint16_t start;
    uint16_t count;

    if (lb->full) {
        start = lb->head;
        count = LAUNCH_BUFFER_SIZE;
    } else {
        start = 0;
        count = lb->head;
    }

    for (uint16_t i = 0; i < count; i++) {
        uint16_t cur_entry = (start + i) % LAUNCH_BUFFER_SIZE;
        write_sd(&data_file, lb->buffer[cur_entry]);
    }

    flush_file(&data_file);

    lb->flushed = 1;
}
