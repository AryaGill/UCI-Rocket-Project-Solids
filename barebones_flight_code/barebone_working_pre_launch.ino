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
#include <Adafruit_AHRS.h>

// #include "air_brakes_drag.h"

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
Adafruit_Mahony algo;

//CSV File Declaration
File dataFile;
File stateFile;
File simFile;

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
unsigned long launch_start_time = 0;
unsigned long drogue_primary_start_time = 0;
unsigned long drogue_primary_end_time = 0;
unsigned long drogue_secondary_start_time = 0;
unsigned long main_primary_start_time = 0;
unsigned long main_primary_end_time = 0;
unsigned long main_secondary_start_time = 0;

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


// Altitude Filtering for Flight State
#define LAUNCH_THRESHOLD 10
#define APOGEE_THRESHOLD -2
#define LANDED_THRESHOLD -0.2
#define ALT_DIF_BUF_SIZE 10
float alt_dif_buffer[ALT_DIF_BUF_SIZE];
int alt_dif_buffer_idx = 0;

//Variables for setting how many data entries are written at once
int cycles_per_write = 20;
int write_count=0;

// Flight State Variables
enum FlightState {
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


void initialize_dataFile() {
  dataFile = SD.open("rocket.csv", FILE_READ);

  if (!dataFile) {
    // File doesn't exist 
    Serial.println("Creating data file.");
    dataFile = SD.open("rocket.csv", FILE_WRITE);
    if (!dataFile) {
        Serial.println("Failed to create file");
    }
  }


  String dataString = "Cam1,Cam2,Temp,Press,Alt,Accel_x2,Accel_y2,Accel_z2,Accel_x,Accel_y,Accel_z,Vel_x,Vel_y,Vel_z,Vel_x2,Vel_y2,Vel_z2,Gyro_x,Gyro_y,Gyro_z,Mag_x,Mag_y,Mag_z,Quaternion_1,Quaternion_2,Quaternion_3,Quaternion_4,Time,State";
  dataFile.println(dataString);
  dataFile.flush();
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
    Accel_y2 = -LIS3DH_data.accel_z;
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
        Accel_y = -LSM9DS1_data.accel_z;
        Accel_z = -LSM9DS1_data.accel_y;

        Gyro_x = -LSM9DS1_data.gyro_x + .0448;
        Gyro_y = -LSM9DS1_data.gyro_z -.0283;
        Gyro_z = -LSM9DS1_data.gyro_y + .0956;

        Mag_x = -LSM9DS1_data.mag_x;
        Mag_y = -LSM9DS1_data.mag_z;
        Mag_z = -LSM9DS1_data.mag_y;

        algo.update(Gyro_x, Gyro_y, Gyro_z, Accel_x, Accel_y, Accel_z, Mag_x, Mag_y, Mag_z);

        float qw, qx, qy, qz;
        algo.getQuaternion(&qw, &qx, &qy, &qz);

        // Store quaternion values
        Quaternion_1 = qw;
        Quaternion_2 = qx;
        Quaternion_3 = qy;
        Quaternion_4 = qz;

      }
  else {
    Serial.println("Failed to get LSM9DS1 data");
  }
  current_time = millis();

  float dt_vel = (current_time - prev_vel_time) / 1000.0f;

  if (dt_vel > 0 && dt_vel < 0.2f) {

      // Integrate LSM9DS1
      Vel_x  += Accel_x  * dt_vel;
      Vel_y  += Accel_y  * dt_vel;
      Vel_z  += Accel_z  * dt_vel;

      // Integrate LIS3DH
      Vel_x2 += Accel_x2 * dt_vel;
      Vel_y2 += Accel_y2 * dt_vel;
      Vel_z2 += Accel_z2 * dt_vel;
  }

  prev_vel_time = current_time;
}

// void init_sim_file(){
//   simFile = SD.open("sim_data.csv", FILE_READ);
 
//   if (!simFile) {
//     Serial.println("Simulation file does not exist");
//     return;
//   }
 
//   // Read header
//   String header = simFile.readStringUntil('\n');
// }

// bool read_sim_data() {
//     delay(10);
//     if (!simFile.available()) {
//         Serial.println("Reached end of file");
//         return false;  // end of sim data
//     }

//     String line = simFile.readStringUntil('\n');
//     line.trim();

//     // If the line itself is "END", stop immediately
//     if (line == "END") {
//         Serial.println("END detected");
//         return false;
//     }

//     // Parse the 5th column (index 4) for Alt
//     int start = 0;
//     int end = -1;
//     for (int i = 0; i <= 4; i++) {
//         start = end + 1;
//         end = line.indexOf(',', start);

//         // If there is no comma (e.g., malformed row), stop safely
//         if (end == -1 && i < 4) {
//             Serial.println("Malformed line, stopping");
//             return false;
//         }
//     }

//     String altStr = line.substring(start, end);
//     altStr.trim();

//     // Check if altitude column says END
//     if (altStr.equalsIgnoreCase("END")) {
//         Serial.println("END detected in ALT column");
//         return false;
//     }

//     Alt = altStr.toFloat();
//     Serial.println(Alt);

//     return true;  // valid data row
// }


float get_avg_alt_dif() {
  float sum = 0;
 
  float largest1 = -1e9, largest2 = -1e9;
  float smallest1 =  1e9, smallest2 =  1e9;
 
  for (int i = 0; i < ALT_DIF_BUF_SIZE; ++i) {
    float v = alt_dif_buffer[i];
    sum += v;
 
    // Update largest values
    if (v > largest1) {
      largest2 = largest1;
      largest1 = v;
    } else if (v > largest2) {
      largest2 = v;
    }
 
    // Update smallest values
    if (v < smallest1) {
      smallest2 = smallest1;
      smallest1 = v;
    } else if (v < smallest2) {
      smallest2 = v;
    }
  }
 
  return (sum - largest1 - largest2 - smallest1 - smallest2) /
         (ALT_DIF_BUF_SIZE - 4);
}

void update_alt_dif_buf(float new_alt_dif) {
  float cur_time = millis();
  if (cur_time == prev_time){
    return;
  }
  alt_dif_buffer[alt_dif_buffer_idx] = new_alt_dif / (cur_time - prev_alt_time) * 1000;
  alt_dif_buffer_idx = (alt_dif_buffer_idx + 1) % ALT_DIF_BUF_SIZE;
  prev_alt_time = cur_time;
}

void set_flight_state(FlightState new_state) {
  flight_state = new_state;

  //Write flight state to state file
  if (SD.exists("rocket_state.csv")){
    SD.remove("rocket_state.csv");
  }

  stateFile = SD.open("rocket_state.csv", FILE_WRITE);
  if (stateFile) {
      stateFile.println(String(startAlt, 8));
      stateFile.flush();
      stateFile.println((int)flight_state);
      stateFile.flush();
  } else {
      Serial.println("Failed to create state file");
  }
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

  // Open rocket state file
  stateFile = SD.open("rocket_state.csv", FILE_READ);

  if (stateFile) {
    if (stateFile.size() > 0) {
        // File exists and has data
        float csv_startAlt = stateFile.readStringUntil('\n').trim().toFloat();
        
        if (startAlt - csv_startAlt > 50){
          flight_state = static_cast<FlightState>(stateFile.readStringUntil('\n').trim().toInt());
          startAlt = csv_startAlt;
          Serial.print("startAlt loaded from file: ");
          Serial.println(startAlt);
          Serial.print("Flight state loaded from file: ");
          Serial.println(state_to_string(flight_state));
        }
        else{
          set_flight_state(LAUNCH_PAD);
        }

    } else {
        // File exists but empty, close previous read mode
        Serial.println("Writing starting altitude.");
        stateFile.close();
        set_flight_state(LAUNCH_PAD);
    }
  } else {
    // File doesn't exist 
    Serial.println("Creating and writing starting altitude.");
    set_flight_state(LAUNCH_PAD);
  }
}

void update_flight_state() {
  update_alt_dif_buf(Alt - pre_alt);

  // Determine Next State
  switch(flight_state) {
    case LAUNCH_PAD:
      // Detect if launched
      if (get_avg_alt_dif() > LAUNCH_THRESHOLD) {
        dataFile.println("LAUNCHED");
        Serial.println("LAUNCHED");

        launch_start_time = millis();
        set_flight_state(MOTOR_BURN);
      }

      break;

    case MOTOR_BURN:
      // Add logic: wait for certain delay or for acceleration to change
      if (millis() - launch_start_time > 1500){
        set_flight_state(GLIDING_ASCENT);
        Serial.println("Motor Burnout");
      }

      break;

    case GLIDING_ASCENT:
      if (get_avg_alt_dif() < APOGEE_THRESHOLD) {
        set_flight_state(DROGUE_PRIMARY_DEPLOYING);
        digitalWrite(drogue_1, HIGH);
        dataFile.println("Primary Drogue Deployed");
        Serial.println("Primary Drogue Deployed");

        
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
        Serial.println("Secondary Drogue Deployed");

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
      if (Alt - startAlt < 229 && Alt - startAlt > 77){
        digitalWrite(main_1, HIGH);
        dataFile.println("Primary Main Deployed");
        Serial.println("Primary Main Deployed");

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
        Serial.println("Secondary Main Deployed");

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

String state_to_string(FlightState state) {
  switch(state) {
    case LAUNCH_PAD:
      return "LAUNCH_PAD";
      break;
    case MOTOR_BURN:
      return "MOTOR_BURN";
      break;
    case GLIDING_ASCENT:
      return "GLIDING_ASCENT";
      break;
    case DROGUE_PRIMARY_DEPLOYING:
      return "DROGUE_PRIMARY_DEPLOYING";
      break;
    case DROGUE_PRIMARY_DEPLOYED:
      return "DROGUE_PRIMARY_DEPLOYED";
      break;
    case DROGUE_SECONDARY_DEPLOYING:
      return "DROGUE_SECONDARY_DEPLOYING";
      break;
    case DROGUE_SECONDARY_DEPLOYED:
      return "DROGUE_SECONDARY_DEPLOYED";
      break;
    case MAIN_PRIMARY_DEPLOYING:
      return "MAIN_PRIMARY_DEPLOYING";
      break;
    case MAIN_PRIMARY_DEPLOYED:
      return "MAIN_PRIMARY_DEPLOYED";
      break;
    case MAIN_SECONDARY_DEPLOYING:
      return "MAIN_SECONDARY_DEPLOYING";
      break;
    case MAIN_SECONDARY_DEPLOYED:
      return "MAIN_SECONDARY_DEPLOYED";
      break;
    case LANDED:
      return "LANDED";
      break;
    default:
      return "UNKNOWN_STATE";
      break;
  }
}

void log_data() {
  int voltage_left = analogRead(camera1_adc);
  int voltage_right = analogRead(camera2_adc);

  String storageDataString = String(voltage_left) + "," + String(voltage_right) + "," + String(Temp, 7) + "," + String(Press, 7) + "," + String(Alt, 7) + "," +
              String(Accel_x2, 7) + "," + String(Accel_y2, 7) + "," + String(Accel_z2, 7) + "," +
              String(Accel_x, 7) + "," + String(Accel_y, 7) + "," + String(Accel_z, 7) + "," +
              String(Vel_x, 7) + "," + String(Vel_y, 7) + "," + String(Vel_z, 7) + "," +
              String(Vel_x2, 7) + "," + String(Vel_y2, 7) + "," + String(Vel_z2, 7) + "," +
              String(Gyro_x, 7) + "," + String(Gyro_y, 7) + "," + String(Gyro_z, 7) + "," +
              String(Mag_x, 7) + "," + String(Mag_y, 7) + "," + String(Mag_z, 7) + "," +
              String(Quaternion_1, 7) + "," + String(Quaternion_2, 7) + "," + 
              String(Quaternion_3, 7) + "," + String(Quaternion_4, 7) + "," +
              String(millis()) + "," + state_to_string(flight_state);

  dataFile.println(storageDataString);
  Serial.println(storageDataString);


  write_count++;

  String dataString = String(voltage_left) + "," + String(voltage_right) + "," + String(Temp, 1) + "," + String(Press, 1) + "," + String(Alt, 1) + "," +
                String(Accel_x2, 1) + "," + String(Accel_y2, 1) + "," + String(Accel_z2, 1) + "," +
                String(Accel_x, 1) + "," + String(Accel_y, 1) + "," + String(Accel_z, 1) + "," +
                String(Gyro_x, 1) + "," + String(Gyro_y, 1) + "," + String(Gyro_z, 1) + "," + String(stage);             

  if(millis() - prev_time > 500){
    HWSERIAL.println(dataString);
    prev_time = millis();
  }


  if(write_count>=cycles_per_write){
    dataFile.flush();
    write_count=0;
  }
}

// void handle_rf_commands() {
//   //RF command handling
//   if(HWSERIAL.available() > 0) {
//     String receivedData = HWSERIAL.readStringUntil('\n');
//     receivedData.trim();
//     Serial.println(receivedData);
//     if(receivedData == "ON"){
//       Serial.println("Camera On Recieved");
//       HWSERIAL.println("TEENSY Camera on");
//       digitalWrite(camera1,HIGH);
//       digitalWrite(camera2,HIGH);
//     }else if (receivedData == "OFF"){
//       Serial.println("Camera Off Recieved");
//       HWSERIAL.println("TEENSY Camera off");
//       digitalWrite(camera1, LOW);
//       digitalWrite(camera2, LOW);
//     }else if (receivedData == "Fire Main P"){
//       Serial.println("Main Primary"); 
//       HWSERIAL.println("TEENSY Fired Main Primary");
//       digitalWrite(main_1, HIGH);
//       delay(charge_delay);
//       digitalWrite(main_1, LOW);
//     } else if (receivedData == "Fire Main S"){
//       Serial.println("Main Secondary"); 
//       HWSERIAL.println("TEENSY Fired Main Secondary");
//       digitalWrite(main_2, HIGH);
//       delay(charge_delay);
//       digitalWrite(main_2, LOW);
//     } else if (receivedData == "Fire Drogue P"){
//       Serial.println("Drogue Primary"); 
//       HWSERIAL.println("TEENSY Fired Drogue Primary");
//       dataFile.println("Drouge Primary");
//       digitalWrite(drogue_1, HIGH);
//       delay(charge_delay);
//       digitalWrite(drogue_1, LOW);
//     } else if (receivedData == "Fire Drogue S"){
//       Serial.println("Drogue Secondary"); 
//       HWSERIAL.println("TEENSY Fired Drogue Secondary");
//       dataFile.println("Drouge Secondary");
//       digitalWrite(drogue_2, HIGH);
//       delay(charge_delay);
//       digitalWrite(drogue_2, LOW); 
//     }  else if (receivedData == "CAM1ON"){
//       Serial.println("Camera1 On Recieved");
//       HWSERIAL.println("TEENSY Camera1 on");
//       digitalWrite(camera1,HIGH);
//     } else if (receivedData == "CAM2ON"){
//       Serial.println("Camera2 On Recieved");
//       HWSERIAL.println("TEENSY Camera2 on");
//       digitalWrite(camera2,HIGH);
//     } else if (receivedData == "CAM1OFF"){
//       Serial.println("Camera1 Off Recieved");
//       HWSERIAL.println("TEENSY Camera1 OFF");
//       digitalWrite(camera1,LOW);
//     } else if (receivedData == "CAM2OFF"){
//       Serial.println("Camera2 Off Recieved");
//       HWSERIAL.println("TEENSY Camera2 OFF");
//       digitalWrite(camera2,LOW);
//     }
//   }
// }

void setup() {
  Serial.begin(115200);

  Serial.println("Running the Flight Computer\n");
  
  pinMode(LED_BUILTIN, OUTPUT);
  digitalWrite(LED_BUILTIN, HIGH);

  HWSERIAL.begin(57600);

  pinMode(main_1, OUTPUT);
  pinMode(main_2, OUTPUT);
  pinMode(drogue_1, OUTPUT);
  pinMode(drogue_2, OUTPUT);
  pinMode(buzzer, OUTPUT);
  pinMode(camera1, OUTPUT);
  pinMode(camera2, OUTPUT);

//Turn both cameras on by default
  digitalWrite(camera1, HIGH);
  digitalWrite(camera2, HIGH);

  initialize_sensors();

  algo.begin(200);

  if (!SD.begin(BUILTIN_SDCARD)) {
    Serial.println("SD card failed or not present.");
  }
  initialize_flight_state();
  initialize_dataFile();

  delay(1000);
  // analogWriteFrequency(buzzer, 4500);
  // analogWrite(buzzer, 128);
  prev_time = millis();

  //comment out for actual launch
  // digitalWrite(buzzer, LOW);

  // initialize_kalman_filter();
  // init_sim_file();
}

void loop(){
  
  
  read_sensors();

  // if (!read_sim_data()) {
  //       Serial.println("Simulation finished. Halting.");
  //       while (true) { delay(100); }  // freeze program
  //   }

  update_flight_state();

  log_data();

  // handle_rf_commands();
}
