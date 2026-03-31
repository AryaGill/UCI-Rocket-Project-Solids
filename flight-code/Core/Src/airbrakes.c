#include "airbrakes.h"
#include "telemetry.h"
#include <math.h>
#include <stdio.h>

// Servo timer
extern TIM_HandleTypeDef htim3;

// Air Brakes variables
#define GAMMA 1.4
#define R 287.05287
#define g 9.80665 // Gravity
#define L 0.0065 // Temperature Lapse Rate
#define DESIRED_SEARCH_TIME 20 // ms
#define TIME_PER_SIM_STEP 0.0132 // ms
float deltaT_coefficient = TIME_PER_SIM_STEP * log2f(NUM_DEPLOYMENT_LEVELS) / DESIRED_SEARCH_TIME / g;
float ground_temp = 300;

float TARGET_APOGEE_M = TARGET_APOGEE_FT * 0.3048;

// Deployment levels should be evenly spread between least and most deployment (inclusive)
// Mach numbers should be evenly spread between 0 and 0.7 (inclusive)
// deployment levels: 0.0, 0.1, 0.2, ..., 1.0
// mach levels: 0.05, 0.1, 0.15, ..., 0.7

// Night Fury CdA
float air_brakes_CdA[NUM_RECORDED_DEPLOYMENT_LEVELS][NUM_RECORDED_MACH_NUMS] = {

	{0.0089884,0.0085976,0.0084022,0.0084022,0.0084022,0.0084022,0.0084022,
	0.0084022,0.0084022,0.0085976,0.0085976,0.008793,0.0089884,0.0091838},

	{0.00943484,0.00904404,0.008853104,0.008853104,0.008857569,0.008866498,0.008870962,
	0.008875426,0.008875426,0.009075291,0.009079755,0.009284084,0.009488413,0.009688277},

	{0.00988128,0.00949048,0.009304009,0.009304009,0.009312938,0.009330795,0.009339724,
	0.009348653,0.009348653,0.009552982,0.00956191,0.009775168,0.009988426,0.010192754},

	{0.0103277,0.0099369,0.009754893,0.009754893,0.009768286,0.009795072,0.009808465,
	0.009821858,0.009821858,0.010030651,0.010044044,0.01026623,0.010488416,0.010697209},

	{0.0107741,0.0103833,0.010205757,0.010205757,0.010223614,0.010259328,0.010277185,
	0.010295042,0.010295042,0.010508299,0.010526156,0.01075727,0.010988384,0.011201641},

	{0.0112206,0.0108298,0.010656722,0.010656722,0.010679044,0.010723688,0.01074601,
	0.010768332,0.010768332,0.010986054,0.011008376,0.01124842,0.011488464,0.011706186},

	{0.011667,0.0112762,0.011107586,0.011107586,0.011134372,0.011187944,0.01121473,
	0.011241516,0.011241516,0.011463702,0.011490488,0.01173946,0.011988432,0.012210618},

	{0.0121135,0.0117227,0.011558551,0.011558551,0.011589802,0.011652304,0.011683555,
	0.011714806,0.011714806,0.011941457,0.011972708,0.01223061,0.012488512,0.012715163},

	{0.0125599,0.0121691,0.012009415,0.012009415,0.01204513,0.01211656,0.012152275,
	0.01218799,0.01218799,0.012419105,0.01245482,0.01272165,0.01298848,0.013219595},

	{0.0130063,0.0126155,0.012460279,0.012460279,0.012500458,0.012580816,0.012620995,
	0.012661174,0.012661174,0.012896753,0.012936932,0.01321269,0.013488448,0.013724027},

	{0.0134528,0.013062,0.012911244,0.012911244,0.012955888,0.013045176,0.01308982,
	0.013134464,0.013134464,0.013374508,0.013419152,0.01370384,0.013988528,0.014228572}

};

static inline float clampf(float x, float min, float max)
{
    if (x < min) return min;
    if (x > max) return max;
    return x;
}

float angle_from_vertical(Telemetry_t *telemetry)
{
    float c = 1.0f - 2.0f * (telemetry->q1*telemetry->q1 + telemetry->q2*telemetry->q2);

    // Clamp for numerical safety
    c = clampf(c, -1.0f, 1.0f);

    return acosf(c);   // radians
}

float get_CdA(uint8_t deployment_level, float mach)
{
	float deployment = ((float)deployment_level) / (float)(NUM_DEPLOYMENT_LEVELS - 1);
    // Clamp inputs to valid range
    deployment = clampf(deployment, 0.0f, 1.0f);
    mach = clampf(mach, 0.05f, 0.7f);

    // Convert to fractional index space
    // deployment step = 0.1
    // mach step = 0.05

    float dep_index_f = deployment * 10.0f;
    float mach_index_f = (mach - 0.05f) * 20.0f;

    uint32_t dep_index = (uint32_t)dep_index_f;
    uint32_t mach_index = (uint32_t)mach_index_f;

    // Prevent overflow at upper edge
    if (dep_index >= NUM_RECORDED_DEPLOYMENT_LEVELS - 1)
        dep_index = NUM_RECORDED_DEPLOYMENT_LEVELS - 2;

    if (mach_index >= NUM_RECORDED_MACH_NUMS - 1)
        mach_index = NUM_RECORDED_MACH_NUMS - 2;

    float dep_frac = dep_index_f - dep_index;
    float mach_frac = mach_index_f - mach_index;

    // Fetch 4 surrounding points
    float CdA00 = air_brakes_CdA[dep_index][mach_index];
    float CdA10 = air_brakes_CdA[dep_index+1][mach_index];
    float CdA01 = air_brakes_CdA[dep_index][mach_index+1];
    float CdA11 = air_brakes_CdA[dep_index+1][mach_index+1];

    // Bilinear interpolation
    float CdA0 = CdA00 + dep_frac * (CdA10 - CdA00);
    float CdA1 = CdA01 + dep_frac * (CdA11 - CdA01);

    float CdA = CdA0 + mach_frac * (CdA1 - CdA0);

    return CdA;
}

float get_mach_number(const float velocity, const float temp){
	float speed_of_sound = pow(R * GAMMA * temp, 0.5);
	return velocity / speed_of_sound;
}

float get_mag2(float x, float y){
	return sqrtf(x*x + y*y);
}

float get_mag3(float x, float y, float z){
	return sqrtf(x*x + y*y + z*z);
}

float predict_apogee(Telemetry_t *telemetry, uint8_t deployment_level){
	// Get time step
	float deltaT = clampf(telemetry->velocity_world_z * deltaT_coefficient, 0.01, 0.1);

	// Get angle from vertical
	float theta = angle_from_vertical(telemetry);
	if(theta > 80.0f * M_PI / 180) theta = 80.0f * M_PI / 180;

	// Initial conditions
	float alt_sim = telemetry->alt_fused;
	float vz_sim = telemetry->velocity_world_z + VELOCITY_TIME_SHIFT * telemetry->accel_world_z;
	float vx_sim = vz_sim * tanf(theta);
//	float vx_sim = get_mag2(telemetry->velocity_world_x, telemetry->velocity_world_y);

	// Convert pressure from hPa to Pa
	float pressure_Pa = telemetry->pressure * 100.0f;
	float temperature_K = fmaxf(ground_temp - (L * telemetry->altitude), 1);

	for (int i = 0; i < 100000; ++i){
		float vz_sim_before = vz_sim;
		float T_local = fmaxf(temperature_K - (L * (alt_sim - telemetry->altitude)), 1);
		float mach_number = get_mach_number(get_mag2(vz_sim, vx_sim), T_local);

		float angle_sim = atan2(vx_sim, vz_sim);
		float airbrake_CdA = angle_sim < 30 * M_PI / 180 ? get_CdA(deployment_level, mach_number) : get_CdA(0, mach_number);

		float p_local = pressure_Pa * pow(T_local / temperature_K, g / (R * L));
		float rho_sim = p_local / (R * T_local);

		float Fd = 0.5 * airbrake_CdA * rho_sim * (vx_sim * vx_sim + vz_sim * vz_sim);

		float Fx = -Fd * sinf(angle_sim);
		float Fz = -Fd * cosf(angle_sim) - g * MASS;
		vx_sim += (Fx / MASS) * deltaT;
		vz_sim += (Fz / MASS) * deltaT;
		alt_sim += ((vz_sim + vz_sim_before) / 2) * deltaT;

		if (vz_sim < 0){
			break;
		}
	}

	return alt_sim;
}

void set_optimal_deployment(FlightState_t flight_state, Telemetry_t *telemetry){
//	float horizontal_speed = get_mag2(telemetry->velocity_world_x, telemetry->velocity_world_y);
//	float angle_of_attack = atan2f(horizontal_speed, telemetry->velocity_world_z);
	float angle_of_attack = angle_from_vertical(telemetry);
	float local_temp = fmaxf(ground_temp - (L * telemetry->altitude), 1);
	if (flight_state != GLIDING_ASCENT
			|| get_mach_number(telemetry->velocity_world_z, local_temp) > 0.7
			|| angle_of_attack > 30 * M_PI / 180){
		set_airbrakes_deployment_level(telemetry, 0);
		return;
	}

	uint8_t low = 0;
	uint8_t high = NUM_DEPLOYMENT_LEVELS - 1;

	int num_sims = log2f(NUM_DEPLOYMENT_LEVELS);

	float pred_apogee;
	int mid;
	for (int i = 0; i < num_sims; ++i){
		mid = ((high + low) / 2) + 1;
		pred_apogee = predict_apogee(telemetry, mid);
		if (pred_apogee >= TARGET_APOGEE_M){
			low = mid;
		}
		else{
			high = mid - 1;
		}
	}

	// Set predicted apogee variable
	if (low == mid) {
		// No need to recompute predicted apogee
		telemetry->predicted_apogee = pred_apogee;
	}
	else {
		// Need to recompute predicted apogee
		telemetry->predicted_apogee = predict_apogee(telemetry, low);
	}

	// Set deployment level
	set_airbrakes_deployment_level(telemetry, low);
}

void init_airbrakes_servo(){
	HAL_TIM_PWM_Start(&htim3, AIRBRAKES_SERVO_1_CHANNEL);
	HAL_TIM_PWM_Start(&htim3, AIRBRAKES_SERVO_2_CHANNEL);
}

void set_airbrakes_initial_temp(Telemetry_t *telemetry){
	ground_temp = telemetry->temperature + 273.15;
}

void set_airbrakes_servo_angle(float angle)
{
	if (angle < 0.0f) angle = 0.0f;
	if (angle > 180.0f) angle = 180.0f;

	uint32_t pulse =
		SERVO_MIN_US +
		(uint32_t)((SERVO_MAX_US - SERVO_MIN_US) * (angle / 180.0f));

	__HAL_TIM_SET_COMPARE(&htim3, AIRBRAKES_SERVO_1_CHANNEL, pulse);
	__HAL_TIM_SET_COMPARE(&htim3, AIRBRAKES_SERVO_2_CHANNEL, pulse);
}

void set_airbrakes_deployment_level(Telemetry_t *telemetry, uint8_t deployment){
	// TODO: change this to be a map from deployment level to servo angle
	// deployment is int from 0 to NUM_DEPLOYMENT_LEVELS - 1
	if (deployment >= NUM_DEPLOYMENT_LEVELS) deployment = NUM_DEPLOYMENT_LEVELS - 1;

	telemetry->airbrake_deployment = deployment;

	// if deployment > 0.851 (max deployment of current night fury), then set to max
	float real_max_percent_of_theoretical_max = 0.851;
	float t = (float)deployment / ((float)(NUM_DEPLOYMENT_LEVELS - 1) * real_max_percent_of_theoretical_max);
	t = clampf(t, 0, 1);
	float angle = SERVO_ANGLE_NOT_EXTENDED + t * (SERVO_ANGLE_EXTENDED - SERVO_ANGLE_NOT_EXTENDED);
	set_airbrakes_servo_angle(angle);
}

void set_airbrakes_deployment_level_on_rail(Telemetry_t *telemetry, uint8_t deployment){
	// TODO: change this to be a map from deployment level to servo angle
	// deployment is int from 0 to NUM_DEPLOYMENT_LEVELS - 1
	if (deployment >= NUM_DEPLOYMENT_LEVELS) deployment = NUM_DEPLOYMENT_LEVELS - 1;

	telemetry->airbrake_deployment = deployment;

	float t = (float)deployment / (float)(NUM_DEPLOYMENT_LEVELS - 1);
	float angle = SERVO_ANGLE_NOT_EXTENDED + t * (SERVO_ANGLE_BEFORE_RAIL - SERVO_ANGLE_NOT_EXTENDED);
	set_airbrakes_servo_angle(angle);
}

void perform_airbrakes_servo_sequence(Telemetry_t *telemetry){
	for (uint8_t i = 0; i < NUM_DEPLOYMENT_LEVELS; ++i){
		set_airbrakes_deployment_level_on_rail(telemetry, i);
		HAL_Delay(10);
	}
	for (uint8_t i = NUM_DEPLOYMENT_LEVELS; i > 0; --i){
		set_airbrakes_deployment_level_on_rail(telemetry, i);
		HAL_Delay(10);
	}

	set_airbrakes_deployment_level_on_rail(telemetry, 0);
	HAL_Delay(500);
	set_airbrakes_deployment_level_on_rail(telemetry, NUM_DEPLOYMENT_LEVELS/4);
	HAL_Delay(500);
	set_airbrakes_deployment_level_on_rail(telemetry, NUM_DEPLOYMENT_LEVELS/2);
	HAL_Delay(500);
	set_airbrakes_deployment_level_on_rail(telemetry, NUM_DEPLOYMENT_LEVELS*3/4);
	HAL_Delay(500);
	set_airbrakes_deployment_level_on_rail(telemetry, NUM_DEPLOYMENT_LEVELS - 1);
	HAL_Delay(1000);
	set_airbrakes_deployment_level_on_rail(telemetry, 0);
	HAL_Delay(500);
	set_airbrakes_deployment_level_on_rail(telemetry, NUM_DEPLOYMENT_LEVELS - 1);
	HAL_Delay(500);
	set_airbrakes_deployment_level_on_rail(telemetry, 0);
}

void set_target_apogee(Telemetry_t *telemetry) {
	TARGET_APOGEE_M = predict_apogee(telemetry, NUM_DEPLOYMENT_LEVELS / 2);
	char msg[32];
	snprintf(msg, sizeof(msg), "Set Target Apogee to %.2f", TARGET_APOGEE_M);
	write_datafile_message(msg);
}
