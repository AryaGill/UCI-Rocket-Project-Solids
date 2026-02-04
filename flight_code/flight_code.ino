#include <ArduinoEigen.h>
#include <ArduinoEigenDense.h>
#include <ArduinoEigenSparse.h>

#include <math.h>

// Include custom libraries

#include <Arduino.h>
#include <SD.h>

#include "BPM390_Module.h"
#include "LIS3DH_Module.h"
#include "LSM9DS1_Module.h"
// #include <MadgwickAHRS.h>
#include <Adafruit_Sensor_Calibration.h>
// #include <Adafruit_AHRS.h>
#include "madgwick.h"

#include "air_brakes_drag.h"

#define main_1 11     // main primary
#define main_2 10    // main secondary     
#define drogue_1 12   // drogue primary
#define drogue_2 9  // drogue secondary
#define buzzer 25    // buzzer
#define camera1 20 // camera1
#define camera2 15 // camera2
#define camera1_adc 21 //camera 1 adc
#define camera2_adc 14 //camera 2 adc

#define HWSERIAL Serial7 // Hardware Serial Needed for RF

// #include "kalman-filter.hpp" //Kalman Filter setup

Adafruit_BMP3XX bmp;
BPM390_Module bmpModule(bmp);
Adafruit_LIS3DH lis = Adafruit_LIS3DH();
LIS3DH_Module LIS3DHModule(lis);
Adafruit_LSM9DS1 lsm = Adafruit_LSM9DS1();
LSM9DS1_Module LSM9DS1Module(lsm);

// Madgwick filter;
// Adafruit_Mahony algo;

//CSV File Declaration
File dataFile;
File stateFile;

// Declare global variables
float Temp = 0;
float Press = 0;
float Alt = 0;
float startAlt = 0;

float Accel_x2 = 0;
float Accel_y2 = 0;
float Accel_z2 = 0;
float Accel_x = 0;
float Accel_y = 0;
float Accel_z = 0;
float Mag_x = 0;
float Mag_y = 0;
float Mag_z = 0;
float Gyro_x = 0;
float Gyro_y = 0;
float Gyro_z = 0;
float Quaternion_1 = 0;
float Quaternion_2 = 0;
float Quaternion_3 = 0;
float Quaternion_4 = 0;

float Vel_x  = 0.0f;
float Vel_y  = 0.0f;
float Vel_z  = 0.0f;
float Vel_x2 = 0.0f;
float Vel_y2 = 0.0f;
float Vel_z2 = 0.0f;

// Declare rocket stage detection variables
const int delay_time = 10;
const int charge_delay = 500; //500
const int backup_delay = 500; //2500
bool launch_flag = 0;

// Time variables
long launch_accel_detected_time = -1;
unsigned long launch_start_time = 0;
unsigned long drogue_primary_start_time = 0;
unsigned long drogue_primary_end_time = 0;
unsigned long drogue_secondary_start_time = 0;
unsigned long main_primary_start_time = 0;
unsigned long main_primary_end_time = 0;
unsigned long main_secondary_start_time = 0;
unsigned long prev_mag_filter_time = 0;

//drogue and main cooldown variables
unsigned long cooldown_time = 10000; //set to how long cooldown should be (10s)

float pre_alt = 0;
int stage;
int base_alt = 500; // Hard-coded base altitude in emergency cases
int counter = 0;
int prev_time;
int prev_alt_time = 0;
unsigned long prev_vel_time = 0;
unsigned long current_time = 0;
unsigned int negative_accel_counter = 0;

// Altitude Filtering for Flight State
#define LAUNCH_ACCEL_THRESHOLD 40
#define RAIL_DELAY_TIME 250
#define LAUNCH_EVAL_PERIOD_TIME 250
#define LAUNCH_THRESHOLD 10
#define APOGEE_THRESHOLD -0.1
#define LANDED_THRESHOLD -0.2
#define ALT_DIF_BUF_SIZE 10
float alt_dif_buffer[ALT_DIF_BUF_SIZE];
int alt_dif_buffer_idx = 0;

//Variables for setting how many data entries are written at once
int cycles_per_write = 20;
int write_count=0;

// Flight State Variables
enum FlightState {
  DISARMED,
  LAUNCH_PAD,
  MOTOR_BURN,
  GLIDING_ASCENT,
  DROGUE_PRIMARY_DEPLOYING,
  DROGUE_PRIMARY_DEPLOYED,
  DROGUE_SECONDARY_DEPLOYING,
  DROGUE_SECONDARY_DEPLOYED,
  MAIN_PRIMARY_DEPLOYING,
  MAIN_PRIMARY_DEPLOYED,
  MAIN_SECONDARY_DEPLOYING,
  MAIN_SECONDARY_DEPLOYED,
  LANDED
};
FlightState flight_state = LAUNCH_PAD;

// Air Brakes variables
#define TARGET_APOGEE_FT 6561
#define TARGET_APOGEE_M TARGET_APOGEE_FT * 0.3048
#define GAMMA 1.4
#define R 287.05287
#define g 9.80665 // Gravity
#define L 0.0065 // Temperature Lapse Rate
#define MASS 2.562797
#define WANTED_AIRBRAKE_ALG_TIME 30 // ms
#define TIME_PER_AIRBRAKE_CALL 0.0125 // ms
float deltaT = 0.01;
float A = pow(0.1016, 2) * M_PI;
float deltaT_coefficient = (TIME_PER_AIRBRAKE_CALL / WANTED_AIRBRAKE_ALG_TIME) / g;
int deployment = 0;

//Complimentary Filter variables
// float vel_baro = 0.0f;
// float vel_imu = 0.0f;
// float vel_baro_averaged = 0.0f; //from 1 second moving average
// float vel_fused = 0.0f;

// float alt_cf = 0.0f;
// float alt_fused = 0.0f;

// float prev_alt_cf = 0.0f;
unsigned long prev_cf_time = 0;

// const float VEL_IMU_W   = 0.99f;

// const float ALT_CF_W  = 0.95f;

// const float BARO_TAU = 1.0f;

// float velocity = 0.0f;   // vertical velocity (m/s)
// float alpha = 0.95f;

float alt_fused = 0.0f;      // fused altitude (m)
float velocity_fused = 0.0f; // fused vertical velocity (m/s)
float prev_baro_alt = 0.0f;  // previous barometer altitude for velocity calculation
 
// Complementary filter weights
float alpha_velocity = 0.99f;  // 99% IMU integrated velocity
float alpha_altitude = 0.95f;  // 95% integrated fused velocity

//Kalman State Variables
// KalmanFilter::MatA Amatrix;
// KalmanFilter::MatQ Q;
// KalmanFilter::MatP P;
// KalmanFilter::MatB B;
// KalmanFilter::MatC C;
// KalmanFilter::MatR Rmatrix;

// KalmanFilter::VecX x;   // [alt, vel, bias]

// KalmanFilter kf(Amatrix, B, C, Q, Rmatrix, P);
// bool kf_initialized = false;
// unsigned long last_kf_time = 0;

// // Kalman filter tuning
// float r_var = 1.8f;     // barometer noise (trust less if high)
// float sigma_a = 0.3f;   // accelerometer noise (m/s^2)
// float q_bias = 1e-5f;   // bias drift

float get_drag_coefficient(const int& deployment_level, const float& mach_number){
  if (mach_number >= 0.7){
    return air_brakes_drag_coefficient[NUM_RECORDED_DEPLOYMENT_LEVELS - 1][NUM_RECORDED_MACH_NUMS - 1];
  }

  float mach_idx = mach_number * (NUM_RECORDED_MACH_NUMS - 1) / 0.7;
  float deployment_idx = (float)deployment_level * (NUM_RECORDED_DEPLOYMENT_LEVELS - 1) / (NUM_DEPLOYMENT_LEVELS - 1);

  // Integer and fractional parts
  int mach_i = (int)mach_idx;                 // lower index
  float mach_frac = mach_idx - (float)mach_i;      // fractional part
  int deployment_i = (int)deployment_idx;
  float deployment_frac = deployment_idx - (float)deployment_i;

  float low_low = air_brakes_drag_coefficient[deployment_i][mach_i];
  float low_high = air_brakes_drag_coefficient[deployment_i][mach_i + 1];

  if (deployment_idx == floor(deployment_idx)){
    // deployment_level is round. No need for upper value
    return low_low + (low_high - low_low) * mach_frac;
  }

  float high_low = air_brakes_drag_coefficient[deployment_i + 1][mach_i];
  float high_high = air_brakes_drag_coefficient[deployment_i + 1][mach_i + 1];

  // --- Bilinear interpolation ---
  float x = mach_frac;
  float y = deployment_frac;

  float result =
      (1 - x) * (1 - y) * low_low +
          x   * (1 - y) * low_high +
      (1 - x) *     y   * high_low +
          x   *     y   * high_high;

  return result;
}

float get_mach_number(const float& velocity, const float& temp){
  float speed_of_sound = pow(R * GAMMA * temp, 0.5);
  return velocity / speed_of_sound;
}

float predict_apogee(const float& alt, const float& temp0, const float& pressure0, const float& angle_of_attack, const float& speed0, const int& deployment_level){
  deltaT = max(0.01, min(speed0 * deltaT_coefficient * cos(angle_of_attack), 0.1));
  
  float alt_sim = alt;
  float vz_sim = speed0 * cos(angle_of_attack);
  float vx_sim = speed0 * sin(angle_of_attack);

  for (int i = 0; i < 100000; ++i){
    float vz_sim_before = vz_sim;
    float T_local = max(temp0 - (L * (alt_sim - alt)), 1);
    float mach_number = get_mach_number(pow(vz_sim * vz_sim + vx_sim * vx_sim, 0.5), T_local);

    float airbrake_Cd = get_drag_coefficient(deployment_level, mach_number);

    float p_local = pressure0 * pow(T_local / temp0, g / (R * L));
    float rho_sim = p_local / (R * T_local);

    float Fd = 0.5 * airbrake_Cd * rho_sim * A * (vx_sim * vx_sim + vz_sim * vz_sim);
    
    float angle_sim = atan2(vx_sim, vz_sim);
    float Fx = -Fd * sin(angle_sim);
    float Fz = -Fd * cos(angle_sim) - g * MASS;
    vx_sim += (Fx / MASS) * deltaT;
    vz_sim += (Fz / MASS) * deltaT;
    alt_sim += ((vz_sim + vz_sim_before) / 2) * deltaT;

    if (vz_sim < 0){
      break;
    }
  }

  // Serial.print("Altitude: ");
  // Serial.println(alt_sim);
  return alt_sim;
}

int optimal_deployment(const float& alt, const float& temp0, const float& pressure0, const float& angle_of_attack, const float& speed0){
  if (angle_of_attack > 30 * M_PI / 180){
    return 0;
  }

  int low = 0;
  int high = NUM_DEPLOYMENT_LEVELS - 1;

  int num_sims = log2(NUM_DEPLOYMENT_LEVELS);

  for (int i = 0; i < num_sims; ++i){
    int mid = (high + low) / 2;
    if (predict_apogee(alt, temp0, pressure0, angle_of_attack, speed0, mid) > TARGET_APOGEE_M){
      low = mid;
    }
    else{
      high = mid;
    }
  }

  return low;
}

void initialize_dataFile() {
  dataFile = SD.open("rocket.csv", FILE_WRITE);

  // if (!dataFile) {
  //   // File doesn't exist 
  //   Serial.println("Creating data file.");
  //   dataFile = SD.open("rocket.csv", FILE_WRITE);
  //   if (!dataFile) {
  //       Serial.println("Failed to create file");
  //   }
  // }

  //data headers
  String dataString = "Cam1,Cam2,Temp,Press,Alt,alt_fused,Accel_x2,Accel_y2,Accel_z2,Accel_x,Accel_y,Accel_z,Vel_x,Vel_y,Vel_z,Vel_x2,Vel_y2,Vel_z2,Gyro_x,Gyro_y,Gyro_z,Mag_x,Mag_y,Mag_z,Quaternion_1,Quaternion_2,Quaternion_3,Quaternion_4,accel_world_x,accel_world_y,accel_world_z,Time,State,Deployment,Predicted_Apogee";
  dataFile.println(dataString);
  dataFile.flush();
}

// void initialize_kalman_filter() {
//   // Initialize Kalman filter
//   float dt = 0.02f;  // initial timestep estimate

//   Amatrix << 1, dt, -0.5f * dt * dt,
//              0, 1,      -dt,
//              0, 0,       1;

//   B << 0.5f * dt * dt,
//        dt,
//        0.0f;

//   // Observation matrix (3 measurements: altitude, vel_y, vel_y2)
//   C << 1, 0, 0,
//        0, 1, 0,
//        0, 1, 0;

//   // Process noise
//   Q.setZero();
//   Q(0,0) = 0.25f * sigma_a * sigma_a * powf(dt,4);
//   Q(1,1) = sigma_a * sigma_a * powf(dt,2);
//   Q(2,2) = q_bias * dt;

//   // Measurement noise
//   Rmatrix.setZero();
//   Rmatrix(0,0) = r_var;  // barometer
//   Rmatrix(1,1) = 2.0f;   // velocity 1
//   Rmatrix(2,2) = 2.0f;   // velocity 2

//   P = KalmanFilter::MatP::Identity() * 100.0f;

//   // Initialize Kalman filter
//   kf = KalmanFilter(Amatrix, B, C, Q, Rmatrix, P);
//   x << startAlt, 0.0f, 0.0f;
//   kf.init(x);
//   kf_initialized = true;

//   last_kf_time = millis();
//   prev_vel_time = millis();
//   current_time = millis();
// }
// Add these global variables at top with other globals
float accel_world_x = 0.0f;  // gravity-compensated vertical acceleration
float accel_world_y = 0.0f;
float accel_world_z = 0.0f;
// Add this function to transform accelerometer to world frame
void transform_accel_to_world() {
  // // Average the two IMU accelerations (body frame)
  // float ax_body = 0.5f * (Accel_x + Accel_x2);
  // float ay_body = 0.5f * (Accel_y + Accel_y2);
  // float az_body = 0.5f * (Accel_z + Accel_z2);
  
  // // Get quaternion from your AHRS filter
  // float qw = Quaternion_1;
  // float qx = Quaternion_2;
  // float qy = Quaternion_3;
  // float qz = Quaternion_4;
  
  // // Rotate acceleration vector from body frame to world frame using quaternion
  // // Formula: v' = q * v * q^(-1)
  // // Simplified for acceleration vector [ax, ay, az]:
  
  // float t2 = qw * qx;
  // float t3 = qw * qy;
  // float t4 = qw * qz;
  // float t5 = -qx * qx;
  // float t6 = qx * qy;
  // float t7 = qx * qz;
  // float t8 = -qy * qy;
  // float t9 = qy * qz;
  // float t10 = -qz * qz;
  
  // float ax_world = 2.0f * ((t8 + t10) * ax_body + (t6 - t4) * ay_body + (t3 + t7) * az_body) + ax_body;
  // float ay_world = 2.0f * ((t4 + t6) * ax_body + (t5 + t10) * ay_body + (t9 - t2) * az_body) + ay_body;
  // float az_world = 2.0f * ((t7 - t3) * ax_body + (t2 + t9) * ay_body + (t5 + t8) * az_body) + az_body;
  
  // // Remove gravity from vertical axis (world frame Z points up, gravity is -9.81 m/s²)
  // accel_world_z = az_world - 9.81f;




  // Average IMUs (body frame)
  float ax = 0.5f * (Accel_x  + Accel_x2);
  float ay = 0.5f * (Accel_y  + Accel_y2);
  float az = 0.5f * (Accel_z  + Accel_z2);

  // Quaternion (w, x, y, z)
  float qw = Quaternion_1;
  float qx = Quaternion_2;
  float qy = Quaternion_3;
  float qz = Quaternion_4;

  // Rotation matrix (body → world)
  float R11 = 1.0f - 2.0f*(qy*qy + qz*qz);
  float R12 = 2.0f*(qx*qy - qz*qw);
  float R13 = 2.0f*(qx*qz + qy*qw);

  float R21 = 2.0f*(qx*qy + qz*qw);
  float R22 = 1.0f - 2.0f*(qx*qx + qz*qz);
  float R23 = 2.0f*(qy*qz - qx*qw);

  float R31 = 2.0f*(qx*qz - qy*qw);
  float R32 = 2.0f*(qy*qz + qx*qw);
  float R33 = 1.0f - 2.0f*(qx*qx + qy*qy);

  // Rotate acceleration into world frame
  float ax_world = R11*ax + R12*ay + R13*az;
  float ay_world = R21*ax + R22*ay + R23*az;
  float az_world = R31*ax + R32*ay + R33*az;

  // Remove gravity (see section below!)
  accel_world_x = ax_world;
  accel_world_y = ay_world;
  accel_world_z = az_world - 9.81f;
}

//complimentary filter
void complementary_filter() {
  unsigned long now = micros();
  float dt = (now - prev_cf_time) * 1e-6f;
  
  if (dt <= 0.0f) return;  // safety check
  
  prev_cf_time = now;
  
  // STAGE 1: Velocity Fusion
  // Use gravity-compensated, tilt-corrected vertical acceleration
  float velocity_imu = velocity_fused + accel_world_z * dt;
  // float velocity_imu = velocity_fused + (-Accel_z - 9.81f) * dt;
  
  // Calculate barometric velocity
  float baro_alt = Alt - startAlt;
  float velocity_baro = (baro_alt - prev_baro_alt) / dt;
  prev_baro_alt = baro_alt;
  
  // Fuse velocities: 99% IMU, 1% barometer
  velocity_fused = alpha_velocity * velocity_imu
                 + (1.0f - alpha_velocity) * velocity_baro;
  
  // STAGE 2: Altitude Fusion
  // Integrate fused velocity to get altitude prediction
  float alt_from_velocity = alt_fused + velocity_fused * dt;
  
  // Fuse altitudes: 95% integrated velocity, 5% raw barometer
  alt_fused = alpha_altitude * alt_from_velocity
            + (1.0f - alpha_altitude) * baro_alt;
}


void initialize_sensors() {
  // BPM390 Setup
  Serial.println("BMP390 Setup");
  if (!bmpModule.begin()){
    Serial.println("Could not find a valid BMP sensor");
    //while (1);
  }

  // LIS3DH Setup
  Serial.println("LIS3DH Setup");
  if (!LIS3DHModule.begin()){
    Serial.println("Could not find a valid LIS3DH sensor");
  }


  // LSM9DS1 Setup
  Serial.println("LSM9DS1 Setup");
  if (!LSM9DS1Module.begin()){
    Serial.println("Could not find a valid LSM9DS1 sensor");
  }

  delay(100);

  // Read accel
  LSM9DS1_SensorData LSM9DS1_data = LSM9DS1Module.readData();
  Accel_x = -LSM9DS1_data.accel_x;
  Accel_y = LSM9DS1_data.accel_z;
  Accel_z = -LSM9DS1_data.accel_y;
}

void read_sensors() {
  // BPM390 Data
  BPM_SensorData BPM_data = bmpModule.readData();
  if (BPM_data.temperature != -999) {
    Temp = BPM_data.temperature;
    Press = BPM_data.pressure;
    Alt = BPM_data.altitude;
  }
  else {
    Serial.println("Failed to get BPM390 data");
  }

  // LIS3DH Data
  LIS3DH_SensorData LIS3DH_data = LIS3DHModule.readData();
  if (LIS3DH_data.accel_x != -999 && LIS3DH_data.accel_y != -999 && LIS3DH_data.accel_z != -999){
    Accel_x2 = LIS3DH_data.accel_x;
    Accel_y2 = LIS3DH_data.accel_z;
    Accel_z2 = -LIS3DH_data.accel_y;
  }
  else {
    Serial.println("Failed to get LIS3DH data");
  }

  // LSM9DS1 Data
  LSM9DS1_SensorData LSM9DS1_data = LSM9DS1Module.readData();

  if (LSM9DS1_data.accel_x != -999 && LSM9DS1_data.accel_y != -999 && LSM9DS1_data.accel_z != -999 &&
      LSM9DS1_data.gyro_x != -999 && LSM9DS1_data.gyro_y != -999 && LSM9DS1_data.gyro_z != -999 &&
      LSM9DS1_data.mag_x != -999 && LSM9DS1_data.mag_y != -999 && LSM9DS1_data.mag_z != -999){
        
        Accel_x = -LSM9DS1_data.accel_x;
        Accel_y = LSM9DS1_data.accel_z;
        Accel_z = -LSM9DS1_data.accel_y;

        Gyro_x = -LSM9DS1_data.gyro_x + .0448;
        Gyro_y = LSM9DS1_data.gyro_z -.0283;
        Gyro_z = -LSM9DS1_data.gyro_y + .0956;

        Mag_x = -LSM9DS1_data.mag_x;
        Mag_y = LSM9DS1_data.mag_z;
        Mag_z = -LSM9DS1_data.mag_y;

        unsigned long cur_time = micros();

        float dt = (cur_time - prev_mag_filter_time) * 1e-6f;
        prev_mag_filter_time = cur_time;

        if (dt > 0){
          Madgwick_Update(Gyro_x, Gyro_y, Gyro_z,
                        Accel_x, Accel_y, Accel_z,
                        Mag_x, Mag_y, Mag_z,
                        dt);
        }

        // algo.update(Gyro_x, Gyro_y, Gyro_z,
        //             Accel_x, Accel_y, Accel_z,
        //             Mag_x, Mag_y, Mag_z,
        //             dt);
        // algo.updateIMU(Gyro_x, Gyro_y, Gyro_z,
        //             Accel_x, Accel_y, Accel_z,
        //             dt);

        // float qw, qx, qy, qz;
        // algo.getQuaternion(&qw, &qx, &qy, &qz);

        // Quaternion_1 = qw;
        // Quaternion_2 = qx;
        // Quaternion_3 = qy;
        // Quaternion_4 = qz;
      }
  else {
    Serial.println("Failed to get LSM9DS1 data");
  }
  current_time = micros();

  float dt_vel = (current_time - prev_vel_time) * 1e-6f;

  if (dt_vel > 0) {

      // Integrate LSM9DS1
      Vel_x  += Accel_x  * dt_vel;
      Vel_y  += Accel_y  * dt_vel;
      Vel_z  += (Accel_z)  * dt_vel;

      // Integrate LIS3DH
      Vel_x2 += Accel_x2 * dt_vel;
      Vel_y2 += Accel_y2 * dt_vel;
      Vel_z2 += (Accel_z2) * dt_vel;
  }

  prev_vel_time = current_time;
}

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
  float cur_time = millis();
  if (cur_time == prev_alt_time){
    return;
  }
  alt_dif_buffer[alt_dif_buffer_idx] = new_alt_dif / (cur_time - prev_alt_time) * 1000;
  alt_dif_buffer_idx = (alt_dif_buffer_idx + 1) % ALT_DIF_BUF_SIZE;
  prev_alt_time = cur_time;
}

void set_flight_state(FlightState new_state) {
  flight_state = new_state;

  // Write flight state to state file
  // stateFile = SD.open("rocket_state.csv", FILE_WRITE);
  // if (stateFile) {
  //     stateFile.println(String(startAlt, 8));
  //     stateFile.flush();
  //     stateFile.println((int)flight_state);
  //     stateFile.flush();
  // } else {
  //     Serial.println("Failed to create state file");
  // }
}

void initialize_flight_state() {
  //possible start altitude after reset fix
  for (int i = 0; i < ALT_DIF_BUF_SIZE + 1; ++i) {
    BPM_SensorData BPM_data = bmpModule.readData();
    if (BPM_data.temperature != -999) {
      startAlt = BPM_data.altitude;
      update_alt_dif_buf(startAlt - pre_alt);
      pre_alt = startAlt;
    }
    else {
    Serial.println("Failed to get BPM390 data");
    }
  }

  set_flight_state(DISARMED);

  // // Open rocket state file
  // stateFile = SD.open("rocket_state.csv", FILE_READ);

  // if (stateFile) {
  //   if (stateFile.size() > 0) {
  //       float read_alt = startAlt;
  //       // File exists and has data
  //       startAlt = stateFile.readStringUntil('\n').trim().toFloat();
  //       Serial.print("startAlt loaded from file: ");
  //       Serial.println(startAlt);

  //       if (read_alt - startAlt > 183){
  //         flight_state = static_cast<FlightState>(stateFile.readStringUntil('\n').trim().toInt());
  //         Serial.print("Flight state loaded from file: ");
  //         Serial.println((int)flight_state);
  //       }
  //       else{
  //         set_flight_state(LAUNCH_PAD);
  //       }

  //   } else {
  //       // File exists but empty, close previous read mode
  //       Serial.println("Writing starting altitude.");
  //       stateFile.close();
  //       set_flight_state(LAUNCH_PAD);
  //   }
  // } else {
  //   // File doesn't exist 
  //   Serial.println("Creating and writing starting altitude.");
  //   set_flight_state(LAUNCH_PAD);
  // }
}

void update_flight_state() {
  // Used for old logic
  update_alt_dif_buf(Alt - pre_alt);

  // Determine Next State
  switch(flight_state) {
    case DISARMED:
      // Rocket Disarmed.
      digitalWrite(LED_BUILTIN, !digitalRead(LED_BUILTIN));
      // analogWrite(buzzer, 0); // Uncomment for testing
      // analogWriteFrequency(buzzer, 4000); // Uncomment for flight
      // analogWrite(buzzer, 128); // Uncomment for flight
      // pinMode(buzzer, OUTPUT);
      break;
    case LAUNCH_PAD:
      // Detect if launched

      // Old logic
      // if (get_avg_alt_dif() > LAUNCH_THRESHOLD) {
        // dataFile.println("LAUNCHED");
        // launch_start_time = millis();
        // set_flight_state(MOTOR_BURN);
      // }

      if (launch_accel_detected_time == -1){
        // Acceleration not detected yet
        //changed to < for vacuum testing
        if (Accel_z > LAUNCH_ACCEL_THRESHOLD){
          // Positive acceleration detected. Begin period of waiting to get off rail.
          launch_accel_detected_time = millis();
          negative_accel_counter = 0;
        }
      }
      else if (millis() - launch_accel_detected_time > RAIL_DELAY_TIME){
        // In evaluation period. Monitor for any negative acceleration value.
        // If detected, reset the system and begin again.
        if (millis() - launch_accel_detected_time > RAIL_DELAY_TIME + LAUNCH_EVAL_PERIOD_TIME){
          // Enough time passed without negative acceleration. Launch detected
          dataFile.println("LAUNCHED");
          launch_start_time = millis();
          set_flight_state(MOTOR_BURN);
        }
        else if (Accel_z < 0){
          // Negative acceleration detected. Reset system.
          if (++negative_accel_counter >= 5){
            launch_accel_detected_time = -1;
          }
        }
      }

      break;

    case MOTOR_BURN:
      // Add logic: wait for certain delay or for acceleration to change
      if (millis() - launch_start_time > 4500){
        set_flight_state(GLIDING_ASCENT);
      }

      break;

    case GLIDING_ASCENT:
      if (get_avg_alt_dif() < APOGEE_THRESHOLD && (Alt - startAlt) > 100) {
        set_flight_state(DROGUE_PRIMARY_DEPLOYING);
        digitalWrite(drogue_1, HIGH);
        dataFile.println("Primary Drogue Deployed");
        //drogue primary starts firing and the time this starts is stored
        drogue_primary_start_time = millis();
      }

      break;

    case DROGUE_PRIMARY_DEPLOYING:
      if (millis() - drogue_primary_start_time >= charge_delay) {
          digitalWrite(drogue_1, LOW);

          //time that primary finishes is stored and bool is set to true so this state does not run again
          drogue_primary_end_time = millis();

          set_flight_state(DROGUE_PRIMARY_DEPLOYED);
        }

      break;

    case DROGUE_PRIMARY_DEPLOYED:
      if (millis() - drogue_primary_end_time >= backup_delay) {
        digitalWrite(drogue_2, HIGH);
        dataFile.println("Secondary Drogue Deployed");

        //time when secondary finishes is stored and bools set so this state does not run again
        drogue_secondary_start_time = millis();

        set_flight_state(DROGUE_SECONDARY_DEPLOYING);
      }

      break;

    case DROGUE_SECONDARY_DEPLOYING:
      if (millis() - drogue_secondary_start_time >= charge_delay){
        digitalWrite(drogue_2, LOW);

        set_flight_state(DROGUE_SECONDARY_DEPLOYED);
      }
    
    case DROGUE_SECONDARY_DEPLOYED:
      // Wait for main deployment
      if (Alt - startAlt < 305 && Alt - startAlt > 77){
        digitalWrite(main_1, HIGH);
        dataFile.println("Primary Main Deployed");
        main_primary_start_time = millis();
        set_flight_state(MAIN_PRIMARY_DEPLOYING);
      }

      break;

    case MAIN_PRIMARY_DEPLOYING:
      if(millis() - main_primary_start_time >= charge_delay){
        digitalWrite(main_1, LOW);
        main_primary_end_time = millis();
        set_flight_state(MAIN_PRIMARY_DEPLOYED);
      }

      break;

    case MAIN_PRIMARY_DEPLOYED:
      if(millis() - main_primary_end_time >= backup_delay){
        digitalWrite(main_2, HIGH);
        dataFile.println("Secondary Main Deployed");
        main_secondary_start_time = millis();
        set_flight_state(MAIN_SECONDARY_DEPLOYING);
      }

      break;

    case MAIN_SECONDARY_DEPLOYING:
      if(millis() - main_secondary_start_time >= charge_delay){
        digitalWrite(main_2, LOW);
        set_flight_state(MAIN_SECONDARY_DEPLOYED);
      }
      break;

    case MAIN_SECONDARY_DEPLOYED:
      if (get_avg_alt_dif() > LANDED_THRESHOLD){ // Change condition
        set_flight_state(LANDED);
      }
      break;

    case LANDED:

      break;
  }

  pre_alt = Alt;
}

// void kalman_filter() {
//   // === Kalman Filter Update ===
//   unsigned long now = millis();
//   float dt = (now - last_kf_time) / 1000.0f;
//   if (dt < 1e-4f) dt = 1e-4f;    // min timestep
//   if (dt > 0.1f)  dt = 0.1f;     // max timestep
//   last_kf_time = now;

//   // update matrices (only A, Q depend on dt)
//   Amatrix << 1, dt, -0.5f*dt*dt,
//       0, 1,      -dt,
//       0, 0,       1;
//   B << 0.5f*dt*dt, dt, 0.0f;
//   Q.setZero();
//   Q(0,0) = 0.25f * sigma_a*sigma_a * powf(dt,4);
//   Q(1,1) = sigma_a*sigma_a * powf(dt,2);
//   Q(2,2) = q_bias * dt;
//   kf.update_dynamics(Amatrix);
//   kf.update_process_noise(Q);

//   // set measurement (baro) and control (accel + gravity)
//   float fused_accel_y = 0.5f * (Accel_y + Accel_y2);  // fuse both Y accelerations
//   Eigen::Matrix<float,1,1> u;
//   u << fused_accel_y + g;

//   // measurement vector (altitude + two Y velocities)
//   Eigen::Matrix<float,3,1> y;
//   y << Alt, Vel_y, Vel_y2;

//   // predict + update
//   kf.predict(u);
//   kf.update(y);
// }

String state_to_string(FlightState state) {
  switch(state) {
    case DISARMED:
      return "0";
      break;
    case LAUNCH_PAD:
      return "1";
      break;
    case MOTOR_BURN:
      return "2";
      break;
    case GLIDING_ASCENT:
      return "3";
      break;
    case DROGUE_PRIMARY_DEPLOYING:
      return "4";
      break;
    case DROGUE_PRIMARY_DEPLOYED:
      return "5";
      break;
    case DROGUE_SECONDARY_DEPLOYING:
      return "6";
      break;
    case DROGUE_SECONDARY_DEPLOYED:
      return "7";
      break;
    case MAIN_PRIMARY_DEPLOYING:
      return "8";
      break;
    case MAIN_PRIMARY_DEPLOYED:
      return "9";
      break;
    case MAIN_SECONDARY_DEPLOYING:
      return "10";
      break;
    case MAIN_SECONDARY_DEPLOYED:
      return "11";
      break;
    case LANDED:
      return "12";
      break;
    default:
      return "13";
      break;
  }
}

void log_data() {
  int voltage_left = analogRead(camera1_adc);
  int voltage_right = analogRead(camera2_adc);

  // retrieve filtered states
  // KalmanFilter::VecX x_hat = kf.state();
  // float Alt_KF = x_hat[0];
  // float Vel_KF = x_hat[1];
  // float Bias_KF = x_hat[2];
// Print combined data

  // String dataString = String(voltage_left) + "," + String(voltage_right) + "," + String(Temp, 7) + "," + String(Press, 7) + "," + String(Alt, 7) + "," +
  //               String(Accel_x2, 7) + "," + String(Accel_y2, 7) + "," + String(Accel_z2, 7) + "," +
  //               String(Accel_x, 7) + "," + String(Accel_y, 7) + "," + String(Accel_z, 7) + "," + String(Alt_KF, 7) + "," + String(Vel_KF, 7) + "," + String(Bias_KF, 7) + "," +
  //               String(Gyro_x, 7) + "," + String(Gyro_y, 7) + "," + String(Gyro_z, 7) + "," +
  //               String(Mag_x, 7) + "," + String(Mag_y, 7) + "," + String(Mag_z, 7) + "," +
  //               String(Quaternion_1, 7) + "," + String(Quaternion_2, 7) + "," + 
  //               String(Quaternion_3, 7) + "," + String(Quaternion_4, 7) + "," +
  //               String(millis()) + "," + state_to_string(flight_state) + "," + String(deployment, 7);

  String storageDataString = String(voltage_left) + "," + String(voltage_right) + "," + String(Temp, 7) + "," + String(Press, 7) + "," + String(Alt-startAlt, 7) + "," + String(alt_fused, 7) + ","+
                String(Accel_x2, 7) + "," + String(Accel_y2, 7) + "," + String(Accel_z2, 7) + "," +
                String(Accel_x, 7) + "," + String(Accel_y, 7) + "," + String(Accel_z, 7) + "," +
                String(Vel_x, 7) + "," + String(Vel_y, 7) + "," + String(Vel_z, 7) + "," +
                String(Vel_x2, 7) + "," + String(Vel_y2, 7) + "," + String(Vel_z2, 7) + "," +
                String(Gyro_x, 7) + "," + String(Gyro_y, 7) + "," + String(Gyro_z, 7) + "," +
                String(Mag_x, 7) + "," + String(Mag_y, 7) + "," + String(Mag_z, 7) + "," +
                String(Quaternion_1, 7) + "," + String(Quaternion_2, 7) + "," + 
                String(Quaternion_3, 7) + "," + String(Quaternion_4, 7) + "," +
                String(accel_world_x, 7) + "," + String(accel_world_y, 7) + "," + String(accel_world_z) + "," +
                String(millis()) + "," + state_to_string(flight_state) + "," + String(deployment) + "," +
                String(predict_apogee(Alt - startAlt, Temp, Press, 0 /*angle of attack*/, 0 /*velocity*/, deployment), 7);


  

  dataFile.println(storageDataString);
  write_count++;



  String dataString = String(millis()) + "," + String(Temp, 1) + "," + String(Press, 1) + "," + String(Alt - startAlt, 1) + "," + String(alt_fused, 1) + "," +
                String(Gyro_x, 1) + "," + String(Gyro_y, 1) + "," + String(Gyro_z, 1) + "," + 
                String(accel_world_x, 1) + "," + String(accel_world_y, 1) + "," + String(accel_world_z, 1) + "," +
                String(Accel_x, 1) + "," + String(Accel_y, 1) + "," + String(Accel_z, 1)+ "," + 
                String(Quaternion_1, 7) + "," + String(Quaternion_2, 7) + "," + 
                String(Quaternion_3, 7) + "," + String(Quaternion_4, 7) + "," +
                state_to_string(flight_state);       

  if(millis() - prev_time > 500){ // CHANGE BACK TO 500
    HWSERIAL.println(dataString);
    Serial.println(dataString);
    prev_time = millis();
  }


  if(write_count>=cycles_per_write){
    dataFile.flush();
    write_count=0;
  }
}

void handle_rf_commands() {
  //RF command handling
  if(HWSERIAL.available() > 0) {
    String receivedData = HWSERIAL.readStringUntil('\n');
    receivedData.trim();
    Serial.println(receivedData);
    if(receivedData == "ON"){
      Serial.println("Camera On Recieved");
      HWSERIAL.println("TEENSY Camera on");
      digitalWrite(camera1,HIGH);
      digitalWrite(camera2,HIGH);

    }else if (receivedData == "OFF"){
      Serial.println("Camera Off Recieved");
      HWSERIAL.println("TEENSY Camera off");
      digitalWrite(camera1, LOW);
      digitalWrite(camera2, LOW);

    }else if (receivedData == "Fire Main P"){
      Serial.println("Main Primary"); 
      HWSERIAL.println("TEENSY Fired Main Primary");
      digitalWrite(main_1, HIGH);
      delay(charge_delay);
      digitalWrite(main_1, LOW);

    } else if (receivedData == "Fire Main S"){
      Serial.println("Main Secondary"); 
      HWSERIAL.println("TEENSY Fired Main Secondary");
      digitalWrite(main_2, HIGH);
      delay(charge_delay);
      digitalWrite(main_2, LOW);

    } else if (receivedData == "Fire Drogue P"){
      Serial.println("Drogue Primary"); 
      HWSERIAL.println("TEENSY Fired Drogue Primary");
      dataFile.println("Drouge Primary");
      digitalWrite(drogue_1, HIGH);
      delay(charge_delay);
      digitalWrite(drogue_1, LOW);

    } else if (receivedData == "Fire Drogue S"){
      Serial.println("Drogue Secondary"); 
      HWSERIAL.println("TEENSY Fired Drogue Secondary");
      dataFile.println("Drogue Secondary");
      digitalWrite(drogue_2, HIGH);
      delay(charge_delay);
      digitalWrite(drogue_2, LOW); 

    }  else if (receivedData == "CAM1ON"){
      Serial.println("Camera1 On Recieved");
      HWSERIAL.println("TEENSY Camera1 on");
      digitalWrite(camera1,HIGH);

    } else if (receivedData == "CAM2ON"){
      Serial.println("Camera2 On Recieved");
      HWSERIAL.println("TEENSY Camera2 on");
      digitalWrite(camera2,HIGH);

    } else if (receivedData == "CAM1OFF"){
      Serial.println("Camera1 Off Recieved");
      HWSERIAL.println("TEENSY Camera1 OFF");
      digitalWrite(camera1,LOW);

    } else if (receivedData == "CAM2OFF"){
      Serial.println("Camera2 Off Recieved");
      HWSERIAL.println("TEENSY Camera2 OFF");
      digitalWrite(camera2,LOW);

    } else if (receivedData == "ARM"){
      Serial.println("Arm Command Recieved");
      HWSERIAL.println("Arm Command Recieved");
      set_flight_state(LAUNCH_PAD);
      // turn off disarmed indicators
      digitalWrite(LED_BUILTIN, HIGH);
      // analogWrite(buzzer, 0); // Uncomment for testing
      analogWriteFrequency(buzzer, 4500); // Uncomment for flight
      analogWrite(buzzer, 128); // Uncomment for flight
    }
  }
}

void setup() {
  Serial.begin(115200);
  // while (!Serial);
  Serial.println("Running the Flight Computer\n");
  
  pinMode(LED_BUILTIN, OUTPUT);
  digitalWrite(LED_BUILTIN, HIGH);
  // pinMode(buzzer, OUTPUT);
  analogWriteFrequency(buzzer, 3500); // Uncomment for flight
  analogWrite(buzzer, 128); // Uncomment for flight

  HWSERIAL.begin(57600);

  pinMode(main_1, OUTPUT);
  pinMode(main_2, OUTPUT);
  pinMode(drogue_1, OUTPUT);
  pinMode(drogue_2, OUTPUT);
  // pinMode(buzzer, OUTPUT);
  pinMode(camera1, OUTPUT);
  pinMode(camera2, OUTPUT);

//Turn both cameras on by default
  digitalWrite(camera1, HIGH);
  digitalWrite(camera2, HIGH);

  initialize_sensors();

  // algo.begin(500);
  Madgwick_Init(&Quaternion_1, &Quaternion_2, &Quaternion_3, &Quaternion_4, Accel_x, Accel_y, Accel_z, 0.1f);

  if (!SD.begin(BUILTIN_SDCARD)) {
    Serial.println("SD card failed or not present.");
  }
  initialize_flight_state();
  initialize_dataFile();

  delay(1000);
  prev_time = millis();
  prev_cf_time = micros();
  prev_mag_filter_time = micros();
  prev_vel_time = micros();

  //comment out for actual launch
  // digitalWrite(buzzer, LOW);
  // analogWriteFrequency(buzzer, 4500); // Uncomment for flight
  // analogWrite(buzzer, 128); // Uncomment for flight

  // initialize_kalman_filter();
}

void loop(){
  read_sensors();
  transform_accel_to_world();
  complementary_filter();

  update_flight_state();

  // kalman_filter();

  // Run Air Brakes Alg
  if (flight_state == GLIDING_ASCENT){ // ADD:  && get_mach_number(velocity, Temp) < 0.7 && angle of attack < 30 deg
    deployment = optimal_deployment(Alt - startAlt, Temp, Press, 0 /*angle of attack*/, 0 /*velocity of roll axis*/);
  }
  else {
    deployment = 0;
  }

  log_data();

  handle_rf_commands();
}
