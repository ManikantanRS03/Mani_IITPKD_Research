// Code to generate a pulse
// Set the set/reset time for pulse by adjusting the duty cycle
// Set the DAC value to set the bias
// Find the switching and detecting bits
// A circular buffer which calculates the probability values

#include "driver/mcpwm.h"
#include "driver/dac.h"
#include "Arduino.h"

#define PWM_PIN 14            // MCPWM output pin
#define PWM_FREQ 50          // PWM frequency in Hz
#define PWM_DUTY 10           // Duty cycle in percentage
#define DAC_PIN DAC_CHANNEL_1 // ESP32 DAC Channel 1 (GPIO 25)
#define ADC_PIN 34            // ADC input pin
#define ADC_THRESHOLD 10    // Example threshold (adjust based on needs)
#define BUFFER_SIZE 500       // Circular buffer size
#define INTERRUPT_PIN 32      // External interrupt pin (connected to PWM)

volatile int adc_value = 0;   // Store ADC reading
volatile int bit_value = 0;   // Store threshold crossing status
volatile int fall_flag = 0;   // Store threshold crossing status
volatile int buffer[BUFFER_SIZE] = {0};  // Circular buffer to store bit values
volatile long int buffer_ktime[BUFFER_SIZE] = {0};  // Circular buffer to store bit values

volatile int buffer_index = 0; // Current index in the buffer
volatile long int buffer_sum = 0;   // Sum of the buffer

volatile int buffer_ktime_index = 0; // Current index in the buffer ktime
volatile long int buffer_ktime_sum = 0;   // Sum of the buffer ktime

volatile float p = 0;
volatile float ktime = 0;

unsigned long start_time = 0;  // Stores time when bit_value = 0
unsigned long end_time = 0;    // Stores time when bit_value = 1
unsigned long switch_time = 0;

volatile int print_count = 0;

// Interrupt Service Routine (ISR) for PWM rising edge
void IRAM_ATTR pwm_rising_isr() {
    // Update circular buffer
    if(bit_value == 1)
    {
      switch_time = end_time - start_time;
      buffer_ktime_sum -= buffer_ktime[buffer_ktime_index]; // Remove old value from ktime buffer
      buffer_ktime[buffer_ktime_index] = switch_time;   // Store new bit value
      buffer_ktime_sum += switch_time;            // Add new value to sum

      buffer_ktime_index = (buffer_ktime_index + 1) % BUFFER_SIZE; // Circular indexing
      ktime = float(buffer_ktime_sum)/BUFFER_SIZE;
      if((buffer_ktime_index == 0) && (print_count < 6))
      {
        if(print_count == 0) {Serial.println("start");}
        Serial.println("p");
        Serial.println(p); // Print sum of buffer
        Serial.println("kt");
        Serial.println(ktime); // Print sum of buffer
        print_count = print_count + 1;
        if(print_count == 6) {Serial.println("end"); print_count = 0;}
      }
    }
    else {
      switch_time = 0;
    }
    buffer_sum -= buffer[buffer_index]; // Remove old value from sum
    buffer[buffer_index] = bit_value;   // Store new bit value
    buffer_sum += bit_value;            // Add new value to sum

    buffer_index = (buffer_index + 1) % BUFFER_SIZE; // Circular indexing
    p = float(buffer_sum)/BUFFER_SIZE;
    bit_value = 0;
    fall_flag = 0;
    return;
}

void setup() {
    Serial.begin(115200);

    // Enable DAC output and set it to ~1V
    dac_output_enable(DAC_PIN);
    dac_output_voltage(DAC_PIN, 255); // 77/255 * 3.3V ≈ 1V

    // Configure MCPWM
    mcpwm_gpio_init(MCPWM_UNIT_0, MCPWM0A, PWM_PIN);
    
    mcpwm_config_t pwm_config;
    pwm_config.frequency = PWM_FREQ;
    pwm_config.cmpr_a = PWM_DUTY;  // Set duty cycle
    pwm_config.cmpr_b = 0;
    pwm_config.counter_mode = MCPWM_UP_COUNTER;
    pwm_config.duty_mode = MCPWM_DUTY_MODE_0;
    mcpwm_init(MCPWM_UNIT_0, MCPWM_TIMER_0, &pwm_config);

    // Configure ADC
    analogReadResolution(12); // Set ADC resolution to 12-bit (0-4095)
    // Attach interrupts
    attachInterrupt(digitalPinToInterrupt(INTERRUPT_PIN), pwm_rising_isr, RISING);
}

void loop() {

    while ((digitalRead(PWM_PIN) == LOW) && (bit_value == 0)){ // Keep checking while in low phase
        if(fall_flag == 0)
        {
          start_time = micros();  // Capture time at falling edge
          fall_flag = 1;
        }
        adc_value = analogRead(ADC_PIN);
        if (adc_value > ADC_THRESHOLD) {
            bit_value = 1;
            end_time = micros();
        }
    }
}
