// Code to generate a pulse
// Set the set/reset time for pulse by adjusting the delay time
// Set the DAC value to set the bias
// Find the switching and detecting bits
// A circular buffer which calculates the probability values

#include "driver/mcpwm.h"
#include "driver/dac.h"
#include "Arduino.h"

#define MOD_PIN 14            // modulation output pin
#define DAC_PIN DAC_CHANNEL_1 // ESP32 DAC Channel 1 (GPIO 25)
#define ADC_PIN 34            // ADC input pin
#define ADC_THRESHOLD 10    // Example threshold (adjust based on needs)
#define BUFFER_SIZE 500       // Circular buffer size

volatile int adc_value = 0;   // Store ADC reading
volatile long int buffer_p[BUFFER_SIZE] = {0};  // Circular buffer to store bit values
volatile int buffer_p_index = 0; // Current index in the buffer ktime
volatile long int buffer_p_sum = 0;   // Sum of the buffer ktime
volatile int bit = 0;
volatile float p = 0;
volatile int print_count = 0;

float ts[19] = {1,2,3,4,6,8,11,15,19,22,27,32,39,46,55,66,79,149,200};
volatile int ts_idx = 0;

void setup() {
    Serial.begin(115200);

    // Enable DAC output and set it to ~1V
    dac_output_enable(DAC_PIN);
    dac_output_voltage(DAC_PIN, 128); // 77/255 * 3.3V ≈ 1V

    // configuring as GPIO
    // Set the GPIO pin as output
    pinMode(MOD_PIN, OUTPUT);
    digitalWrite(MOD_PIN, HIGH);

    // Configure ADC
    analogReadResolution(12); // Set ADC resolution to 12-bit (0-4095)
}

void loop() {
    //while the loop starts we set it as HIGH (reset state)
    if(ts_idx < 20)
    {
    digitalWrite(MOD_PIN, HIGH);
    delay(10); // ms hold in high
    digitalWrite(MOD_PIN, LOW);
    delay(ts[ts_idx]);
    adc_value = analogRead(ADC_PIN);
      if (adc_value > ADC_THRESHOLD) { bit = 1;}
      else {bit = 0;}

      buffer_p_sum -= buffer_p[buffer_p_index]; // Remove old value from p buffer
      buffer_p[buffer_p_index] = bit;   // Store new bit value
      buffer_p_sum += bit;            // Add new value to sum

      buffer_p_index = (buffer_p_index + 1) % BUFFER_SIZE; // Circular indexing
      p = float(buffer_p_sum)/BUFFER_SIZE;

    if(buffer_p_index == 0)
      {
        Serial.println("ts");
        Serial.println(ts[ts_idx]);
        Serial.println("p");
        Serial.println(p); // Print sum of buffer
        print_count = print_count + 1;
        if(print_count == 1) {print_count = 0; ts_idx += 1; if(ts_idx == 19) {Serial.println("end");}}
      }
    }
}
