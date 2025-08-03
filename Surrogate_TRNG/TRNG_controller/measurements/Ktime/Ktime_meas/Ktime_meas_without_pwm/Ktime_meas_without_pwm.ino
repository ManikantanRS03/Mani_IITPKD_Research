// Code to generate a pulse
// Set the set/reset time for pulse by adjusting the duty cycle
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
volatile long int buffer_ktime[BUFFER_SIZE] = {0};  // Circular buffer to store bit values
volatile int buffer_ktime_index = 0; // Current index in the buffer ktime
volatile long int buffer_ktime_sum = 0;   // Sum of the buffer ktime
volatile float ktime = 0;
unsigned long start_time = 0;  // Stores time when bit_value = 0
unsigned long end_time = 0;    // Stores time when bit_value = 1
unsigned long switch_time = 0;
volatile int print_count = 0;


void setup() {
    Serial.begin(115200);

    // Enable DAC output and set it to ~1V
    dac_output_enable(DAC_PIN);
    dac_output_voltage(DAC_PIN, 191); // 77/255 * 3.3V ≈ 1V

    // configuring as GPIO
    // Set the GPIO pin as output
    pinMode(MOD_PIN, OUTPUT);
    digitalWrite(MOD_PIN, HIGH);

    // Configure ADC
    analogReadResolution(12); // Set ADC resolution to 12-bit (0-4095)
}

void loop() {
    //while the loop starts we set it as HIGH (reset state)
   digitalWrite(MOD_PIN, HIGH);
   delay(10); // ms hold in high
   digitalWrite(MOD_PIN, LOW);
   start_time = micros();
   while (true)
   {
    adc_value = analogRead(ADC_PIN);
    if (adc_value > ADC_THRESHOLD)
    {
      end_time = micros();
      switch_time = end_time - start_time;
      buffer_ktime_sum -= buffer_ktime[buffer_ktime_index]; // Remove old value from ktime buffer
      buffer_ktime[buffer_ktime_index] = switch_time;   // Store new bit value
      buffer_ktime_sum += switch_time;            // Add new value to sum

      buffer_ktime_index = (buffer_ktime_index + 1) % BUFFER_SIZE; // Circular indexing
      ktime = float(buffer_ktime_sum)/BUFFER_SIZE;
      break;
    }
   }
   if((buffer_ktime_index == 0) && (print_count < 6))
    {
      if(print_count == 0) {Serial.println("start");}
      Serial.println("kt");
      Serial.println(ktime); // Print sum of buffer
      print_count = print_count + 1;
      if(print_count == 2) {Serial.println("end"); print_count = 0;}
    }
   delay(5); // ms hold in high
}
