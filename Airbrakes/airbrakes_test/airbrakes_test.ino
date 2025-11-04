#include <math.h>

#define TARGET_APOGEE_FT 10000
#define TARGET_APOGEE_M TARGET_APOGEE_FT * 0.3048

float get_drag_coefficient(const int& deployment_level, const float& mach_number){
  return 0.5;
}

float deltaT = 0.01;
float GAMMA = 1.4;
float R = 287.05287;
float g = 9.80665;
float L = 0.0065;
float A = pow(0.07886715773, 2) * M_PI;
float mass = 20;

float wanted_time = 30; // ms
float time_per_call = 0.0125; // ms
float deltaT_coefficient = (time_per_call / wanted_time) / g;
float predict_apogee(const float& alt, const float& temp0, const float& pressure0, const float& angle_of_attack, const float& speed0, const int& deployment_level){
  deltaT = max(0.01, min(speed0 * deltaT_coefficient * cos(angle_of_attack), 0.1));
  Serial.print("DeltaT: ");
  Serial.println(deltaT, 4);
  
  float alt_sim = alt;
  float vz_sim = speed0 * cos(angle_of_attack);
  float vx_sim = speed0 * sin(angle_of_attack);

  for (int i = 0; i < 100000; ++i){
    float vz_sim_before = vz_sim;
    float T_local = max(temp0 - (L * (alt_sim - alt)), 1);
    float speed_of_sound = pow(R * GAMMA * T_local, 0.5);
    float mach_number = pow(vz_sim * vz_sim + vx_sim * vx_sim, 0.5) / speed_of_sound;

    float airbrake_Cd = get_drag_coefficient(deployment_level, mach_number);

    float p_local = pressure0 * pow(T_local / temp0, g / (R * L));
    float rho_sim = p_local / (R * T_local);

    float Fd = 0.5 * airbrake_Cd * rho_sim * A * (vx_sim * vx_sim + vz_sim * vz_sim);
    
    float angle_sim = atan2(vx_sim, vz_sim);
    float Fx = -Fd * sin(angle_sim);
    float Fz = -Fd * cos(angle_sim) - g * mass;
    vx_sim += (Fx / mass) * deltaT;
    vz_sim += (Fz / mass) * deltaT;
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

  int num_sims = 10;

  int low = 0;
  int high = pow(2, num_sims);

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

int deployment = 0;

void setup() {
  pinMode(LED_BUILTIN, OUTPUT);
}

void loop() {
  float alt = 0;
  float temp = 300;
  float pressure = 101325;
  float angle_of_attack = 0;
  float speed = 300;

  unsigned long startTime = millis();
  deployment = optimal_deployment(alt, temp, pressure, angle_of_attack, speed);
  unsigned long endTime = millis();

  Serial.print("Function execution time: ");
  Serial.print(endTime - startTime);
  Serial.println(" milliseconds");

  Serial.print("Deployment Level: ");
  Serial.println(deployment);

  digitalWrite(LED_BUILTIN, HIGH);
  delay(1000);
  digitalWrite(LED_BUILTIN, LOW);
  delay(1000);
}
