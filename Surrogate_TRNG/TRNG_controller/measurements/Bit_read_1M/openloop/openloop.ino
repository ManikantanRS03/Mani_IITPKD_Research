// Code to generate a pulse (PWM), the set reset time is fixed
// Set the DAC value to set the bias (this is the control parameter)
// Find the switching and detecting bits
// A circular buffer which calculates the probability values  (to estimate the probability)
// A PI controller which determine the bias from the error

#include "driver/mcpwm.h"
#include "driver/dac.h"
#include "Arduino.h"

#define PWM_PIN 14            // MCPWM output pin
#define PWM_FREQ 50          // PWM frequency in Hz
#define PWM_DUTY 10           // Duty cycle in percentage
#define INTERRUPT_PIN 32      // External interrupt pin (connected to PWM)

#define DAC_PIN DAC_CHANNEL_1 // ESP32 DAC Channel 1 (GPIO 25)
#define ADC_PIN 34            // ADC input pin
#define ADC_THRESHOLD 10    // Example threshold (adjust based on needs)
#define BUFFER_SIZE 16       // Circular buffer size

volatile int adc_value = 0;   // Store ADC reading
volatile int fall_flag = 0;
volatile long int buffer_p[BUFFER_SIZE] = {0};  // Circular buffer to store bit values
volatile int buffer_p_index = 0; // Current index in the buffer ktime
volatile long int buffer_p_sum = 0;   // Sum of the buffer ktime
volatile int bit = 0;
volatile float p = 0;
volatile long int print_count = 0;


volatile float Ts = 0.02; // sampling time
volatile float set_p = 0.7;
volatile float meas_p = 0;
volatile float error = 0;
volatile float error_prev = 0;

volatile long int count = 0; // for giving a probability set value
volatile int pset_count = 0;

volatile int DAC_val = 80;


void printBufferAsHex() {
    Serial.println("Hex:");
    for (int i = 0; i < BUFFER_SIZE; i += 4) {
        int hexValue = (buffer_p[i] & 1) << 3 |  // MSB
                       (buffer_p[i + 1] & 1) << 2 |
                       (buffer_p[i + 2] & 1) << 1 |
                       (buffer_p[i + 3] & 1);    // LSB
        Serial.printf("%X", hexValue);  // Print as a single hex digit
    }
    Serial.printf("\n");
}

void P_estimate()
{  
    buffer_p_sum -= buffer_p[buffer_p_index]; // Remove old value from sum
    buffer_p[buffer_p_index] = bit;   // Store new bit value
    buffer_p_sum += bit;            // Add new value to sum

    buffer_p_index = (buffer_p_index + 1) % BUFFER_SIZE; // Circular indexing
    // Serial.println(buffer_sum); // Print sum of buffer
    p = float(buffer_p_sum)/BUFFER_SIZE;
    // Serial.println("p");
    // Serial.println(p); // Print sum of buffer
    // Serial.println(switch_time); 
    bit = 0;
    fall_flag = 0;
    return;
}


void IRAM_ATTR pwm_rising_isr() {
  P_estimate();
  count  = count + 1;
  if(count == BUFFER_SIZE)
  {
    printBufferAsHex();
    // Serial.println(p);
    print_count = print_count + 1;
    if(print_count == 125000){Serial.println("end"); print_count = 0;}
    count = 0;
  }
}

void setup() {
    Serial.begin(115200);

    // Enable DAC output and set it to ~1V
    dac_output_enable(DAC_PIN);
    dac_output_voltage(DAC_PIN, DAC_val); // 77/255 * 3.3V ≈ 1V

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
  while ((digitalRead(PWM_PIN) == LOW) && (bit == 0))
  { 
    // Keep checking while in low phase
    if(fall_flag == 0) {fall_flag = 1;}
    adc_value = analogRead(ADC_PIN);
    if (adc_value > ADC_THRESHOLD) {bit = 1;}
  }
}
