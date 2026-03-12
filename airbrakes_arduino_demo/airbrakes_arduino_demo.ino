#include <Servo.h>

#define AIRBRAKES_SERVO_PIN 9

#define SERVO_MIN_US 1000
#define SERVO_MAX_US 2000

#define NUM_DEPLOYMENT_LEVELS 64

#define SERVO_ANGLE_NOT_EXTENDED 180
#define SERVO_ANGLE_EXTENDED 137

Servo airbrakeServo;

void set_airbrakes_servo_pulse(uint16_t pulse)
{
    if (pulse < SERVO_MIN_US) pulse = SERVO_MIN_US;
    if (pulse > SERVO_MAX_US) pulse = SERVO_MAX_US;

    airbrakeServo.writeMicroseconds(pulse);
}

void set_airbrakes_servo_angle(float angle)
{
    if (angle < 0.0f) angle = 0.0f;
    if (angle > 180.0f) angle = 180.0f;

    uint16_t pulse =
        SERVO_MIN_US +
        (uint16_t)((SERVO_MAX_US - SERVO_MIN_US) * (angle / 180.0f));

    set_airbrakes_servo_pulse(pulse);
}

void set_airbrakes_deployment_level(uint8_t deployment)
{
    if (deployment >= NUM_DEPLOYMENT_LEVELS)
        deployment = NUM_DEPLOYMENT_LEVELS - 1;

    float t = (float)deployment / (float)(NUM_DEPLOYMENT_LEVELS - 1);

    float angle =
        SERVO_ANGLE_NOT_EXTENDED +
        t * (SERVO_ANGLE_EXTENDED - SERVO_ANGLE_NOT_EXTENDED);

    set_airbrakes_servo_angle(angle);
}

void perform_airbrakes_servo_sequence()
{
    for (uint8_t i = 0; i < NUM_DEPLOYMENT_LEVELS; ++i)
    {
        set_airbrakes_deployment_level(i);
        delay(10);
    }

    for (int i = NUM_DEPLOYMENT_LEVELS - 1; i >= 0; --i)
    {
        set_airbrakes_deployment_level(i);
        delay(10);
    }

    set_airbrakes_deployment_level(0);
    delay(500);

    set_airbrakes_deployment_level(NUM_DEPLOYMENT_LEVELS/4);
    delay(500);

    set_airbrakes_deployment_level(NUM_DEPLOYMENT_LEVELS/2);
    delay(500);

    set_airbrakes_deployment_level(NUM_DEPLOYMENT_LEVELS*3/4);
    delay(500);

    set_airbrakes_deployment_level(NUM_DEPLOYMENT_LEVELS - 1);
    delay(1000);

    set_airbrakes_deployment_level(0);
    delay(500);

    set_airbrakes_deployment_level(NUM_DEPLOYMENT_LEVELS - 1);
    delay(500);

    set_airbrakes_deployment_level(0);
}

void setup() {
  airbrakeServo.attach(AIRBRAKES_SERVO_PIN);
}

void loop() {
  perform_airbrakes_servo_sequence();
  delay(5000);
}