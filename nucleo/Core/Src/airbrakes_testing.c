#include "airbrakes_testing.h"
#include <math.h>

// Servo timer
extern TIM_HandleTypeDef htim3;

// Air Brakes variables
#define TARGET_APOGEE_FT 8400
#define TARGET_APOGEE_M TARGET_APOGEE_FT * 0.3048
#define GAMMA 1.4
#define R 287.05287
#define g 9.80665 // Gravity
#define L 0.0065 // Temperature Lapse Rate
#define MASS 25.71 // kg
#define DESIRED_SEARCH_TIME 200 // ms
#define TIME_PER_SIM_STEP 0.23 // ms
float deltaT_coefficient = TIME_PER_SIM_STEP * log2f(NUM_DEPLOYMENT_LEVELS) / DESIRED_SEARCH_TIME / g;

// Deployment levels should be evenly spread between least and most deployment (inclusive)
// Mach numbers should be evenly spread between 0 and 0.7 (inclusive)
// deployment levels: 0.0, 0.1, 0.2, ..., 1.0
// mach levels: 0.05, 0.1, 0.15, ..., 0.7

float air_brakes_CdA[NUM_RECORDED_DEPLOYMENT_LEVELS][NUM_RECORDED_MACH_NUMS] =
{
    // 0% (0.02439 m²)
    {0.0070731f,0.00721944f,0.007317f,0.00741456f,0.00743895f,0.00743895f,0.00746334f,0.00748773f,0.00748773f,0.00748773f,0.00751212f,0.00753651f,0.0075609f,0.0075609f},

    // 10% (0.02489 m²)
    {0.007467f,0.00761634f,0.0077159f,0.00776568f,0.00786524f,0.00789013f,0.00793991f,0.0079648f,0.00798969f,0.00801458f,0.00803947f,0.00806436f,0.00806436f,0.00808925f},

    // 20% (0.02539 m²)
    {0.00765355f,0.00786328f,0.00792168f,0.00794961f,0.00797754f,0.00800547f,0.0080334f,0.00806133f,0.00808925f,0.00811972f,0.00814765f,0.00817558f,0.00820351f,0.00823144f},

    // 30% (0.02590 m²)
    {0.00800113f,0.00823566f,0.0083261f,0.00839045f,0.00841109f,0.00843173f,0.00845238f,0.00847304f,0.00849369f,0.00851433f,0.00853497f,0.00855561f,0.00857626f,0.0085969f},

    // 40% (0.02640 m²)
    {0.0082368f,0.008316f,0.0083952f,0.0084744f,0.0085536f,0.0086328f,0.008712f,0.0087912f,0.0088704f,0.0089496f,0.0090288f,0.009108f,0.0091872f,0.0092664f},

    // 50% (0.02690 m²)
    {0.008608f,0.0087963f,0.0089039f,0.0089577f,0.0090115f,0.0090384f,0.0090922f,0.0090922f,0.0090922f,0.0090922f,0.0090922f,0.0090922f,0.0090922f,0.0091191f},

    // 60% (0.02740 m²)
    {0.0088228f,0.00887486f,0.00892618f,0.0089775f,0.00902882f,0.00908014f,0.00913146f,0.00918278f,0.0092341f,0.00928542f,0.00933674f,0.00938806f,0.00943938f,0.0094907f},

    // 70% (0.02790 m²)
    {0.0092907f,0.0095418f,0.0096534f,0.0097092f,0.009765f,0.0098208f,0.0098766f,0.0099324f,0.0099882f,0.010044f,0.0100719f,0.0100998f,0.0101277f,0.0101556f},

    // 80% (0.02841 m²)
    {0.00940371f,0.00968781f,0.00980145f,0.00985827f,0.00991509f,0.00997191f,0.01002873f,0.01008555f,0.01014237f,0.01019919f,0.0102276f,0.01025601f,0.01028442f,0.01031283f},

    // 90% (0.02891 m²)
    {0.0098294f,0.01003177f,0.01023314f,0.01034878f,0.0104066f,0.01046442f,0.01052224f,0.01055115f,0.01058006f,0.01060897f,0.0106957f,0.01078243f,0.01078243f,0.01078243f},

    // 100% (0.02941 m²)
    {0.0102935f,0.01055719f,0.01067583f,0.01076406f,0.01082288f,0.01070524f,0.01091011f,0.01093952f,0.01093952f,0.01096893f,0.01105716f,0.01108657f,0.01111598f,0.01114539f}
};

static inline float clampf(float x, float min, float max)
{
    if (x < min) return min;
    if (x > max) return max;
    return x;
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
	float deltaT = clampf(telemetry->velocity_world_z * deltaT_coefficient, 0.01, 0.1);

	float alt_sim = telemetry->altitude;
//	float vz_sim = telemetry->velocity_r * cosf(telemetry->angle_of_attack);
//	float vx_sim = telemetry->velocity_r * sinf(telemetry->angle_of_attack);
	float vz_sim = telemetry->velocity_world_z;
	float vx_sim = get_mag2(telemetry->velocity_world_x, telemetry->velocity_world_y);

	// Convert pressure from hPa to Pa
	float pressure_Pa = telemetry->pressure * 100.0f;
	float temperature_K = telemetry->temperature + 273.15f;

	for (int i = 0; i < 100000; ++i){
		float vz_sim_before = vz_sim;
		float T_local = fmaxf(temperature_K - (L * (alt_sim - telemetry->altitude)), 1);
		float mach_number = get_mach_number(pow(vz_sim * vz_sim + vx_sim * vx_sim, 0.5), T_local);

		float airbrake_CdA = get_CdA(deployment_level, mach_number);

		float p_local = pressure_Pa * pow(T_local / temperature_K, g / (R * L));
		float rho_sim = p_local / (R * T_local);

		float Fd = 0.5 * airbrake_CdA * rho_sim * (vx_sim * vx_sim + vz_sim * vz_sim);

		float angle_sim = atan2(vx_sim, vz_sim);
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
	float horizontal_speed = get_mag2(telemetry->velocity_world_x, telemetry->velocity_world_y);
	float angle_of_attack = atan2f(horizontal_speed, telemetry->velocity_world_z);
	if (flight_state != GLIDING_ASCENT
			|| get_mach_number(get_mag3(telemetry->velocity_world_x, telemetry->velocity_world_y, telemetry->velocity_world_z), telemetry->temperature + 273.15) > 0.7
			|| angle_of_attack > 30 * M_PI / 180){
		set_airbrakes_deployment_level(telemetry, 0);
		return;
	}

	uint8_t low = 0;
	uint8_t high = NUM_DEPLOYMENT_LEVELS - 1;

	int num_sims = log2f(NUM_DEPLOYMENT_LEVELS);

	float pred_apogee;
	for (int i = 0; i < num_sims; ++i){
		int mid = ((high + low) / 2) + 1;
		pred_apogee = predict_apogee(telemetry, mid);
		if (pred_apogee >= TARGET_APOGEE_M){
			low = mid;
		}
		else{
			high = mid - 1;
		}
	}

	// set predicted apogee variable (wastes time. Comment out if don't want data logged)
	telemetry->predicted_apogee = predict_apogee(telemetry, low);

	// Set deployment level
	set_airbrakes_deployment_level(telemetry, low);
}

void init_airbrakes_servo(){
	HAL_TIM_PWM_Start(&htim3, AIRBRAKES_SERVO_1_CHANNEL);
	HAL_TIM_PWM_Start(&htim3, AIRBRAKES_SERVO_2_CHANNEL);
}

void set_airbrakes_servo_angle(uint8_t angle)
{
    if (angle > 180) angle = 180;

    uint32_t pulse =
        SERVO_MIN_US +
        ((SERVO_MAX_US - SERVO_MIN_US) * angle) / 180;

    __HAL_TIM_SET_COMPARE(&htim3, AIRBRAKES_SERVO_1_CHANNEL, pulse);
    __HAL_TIM_SET_COMPARE(&htim3, AIRBRAKES_SERVO_2_CHANNEL, pulse);
}

void set_airbrakes_deployment_level(Telemetry_t *telemetry, uint8_t deployment){
	// TODO: change this to be a map from deployment level to servo angle
	// deployment is int from 0 to NUM_DEPLOYMENT_LEVELS - 1
	telemetry->airbrake_deployment = deployment;
	set_airbrakes_servo_angle((uint8_t)(((uint32_t)(((float)deployment) / (float)(NUM_DEPLOYMENT_LEVELS - 1) * 180)) / NUM_DEPLOYMENT_LEVELS));
}
