#include <SPI.h>
#include <LoRa.h>
 
#define SS   8
#define RST  4
#define DIO0 7
 
String serialBuffer = "";
 
void setup()
{
    Serial.begin(115200);
 
    while(!Serial);   // wait for USB serial
 
    LoRa.setPins(SS, RST, DIO0);
 
    if (!LoRa.begin(433E6))
    {
        Serial.println("LoRa init failed");
        while (1);
    }
 
    Serial.println("LoRa bridge ready");
}
 
void loop()
{
    readSerialCommands();
    readLoRaTelemetry();
}
 
void readSerialCommands()
{
    while (Serial.available())
    {
        char c = Serial.read();
 
        if (c == '\n')
        {
            sendCommand(serialBuffer);
            serialBuffer = "";
        }
        else
        {
            serialBuffer += c;
        }
    }
}
 
void sendCommand(String cmd)
{
    cmd.trim();
 
    if (cmd.length() == 0)
        return;
 
    LoRa.beginPacket();
    LoRa.print(cmd);
    LoRa.endPacket();
 
    Serial.println(cmd);
}
 
void readLoRaTelemetry()
{
    int packetSize = LoRa.parsePacket();
 
    if (packetSize)
    {
        while (LoRa.available())
        {
            char c = (char)LoRa.read();
            Serial.print(c);
        }
 
        Serial.println();
    }
}