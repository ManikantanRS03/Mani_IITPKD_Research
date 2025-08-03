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
#define BUFFER_SIZE 100       // Circular buffer size

volatile int adc_value = 0;   // Store ADC reading
volatile int fall_flag = 0;
volatile long int buffer_p[BUFFER_SIZE] = {0};  // Circular buffer to store bit values
volatile int buffer_p_index = 0; // Current index in the buffer ktime
volatile long int buffer_p_sum = 0;   // Sum of the buffer ktime
volatile int bit = 0;
volatile float p = 0;
volatile int print_count = 0;


volatile float Ts = 0.02; // sampling time
volatile float set_p = 0.5;
volatile float meas_p = 0;
volatile float error = 0;
volatile float error_prev = 0;

volatile float Kp = 0.2;
volatile float Ki = 0.08;
volatile float integ_bias = 0;
volatile float prop_bias = 0;
volatile float bias = 0;
volatile int int_enb = 1; // intgral windup enable control

volatile long int count = 0; // for giving a probability set value
volatile int pset_count = 0;

volatile int DAC_val = 255;

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

void PI_control()
{
  prop_bias = Kp*error;

  if(int_enb == 1)
  {
    integ_bias = integ_bias + Ki*(error+error_prev)*Ts/2.0;
  }

  // integral windup
  bias = prop_bias + integ_bias;
  //clipping the control output
  if(bias > 0.3)
  {
    bias = 0.3;
    int_enb = 0;
  }
  
  else if(bias < 0 )
  {
    bias = 0;
    int_enb = 0;
  }
  else
  {
    int_enb = 1;
  }

  DAC_val = int(255*(bias*10/3.3));
  error_prev = error;
  if(print_count == 10) {Serial.println(p); Serial.println(bias); print_count = 0;}
  // Serial.println(error);
  // Serial.println(bias);
  // Serial.println(DAC_val);
  dac_output_voltage(DAC_PIN, DAC_val); // 77/255 * 3.3V ≈ 1V
}


void IRAM_ATTR pwm_rising_isr() {
  count  = count + 1;
  print_count = print_count + 1;
  if(count == 2000)
  {
    set_p = 0.5;
    pset_count += 1;
    if(pset_count == 10){Serial.println("end"); pset_count = 0;}
    count = 0;
    Serial.println("set p = ");
    Serial.println(set_p);
  }
  P_estimate();
  meas_p = p;
  error = set_p - meas_p;
  PI_control();

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
