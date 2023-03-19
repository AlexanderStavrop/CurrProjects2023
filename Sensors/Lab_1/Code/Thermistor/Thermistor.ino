#include <math.h>             // Math library for log function
// Thermistor constants
#define beta 4250             // Thermistor beta value 
#define T_0 298.0f            // Thermistor temperature at 25C
// Circuit constants
#define v_in 5.0f             // Arduino input voltage
#define R_2 100000.0f         // Resistor 2 value
// I/O constants
#define ANALOG_PIN A0         // Thermistor analog pin
#define LED_PIN 8             // LED digital pin
// Experiment constants
#define KELVIN_ZERO 273.2f    // Kelvin zero
#define TARGET_TEMP 30.0f     // Target temperature
#define DELAY 500             // Delay between measurements


//.......................................Setup........................................
void setup() {
  Serial.begin(9600);         // Starting the serial communication  
  pinMode(LED_PIN, OUTPUT);   // Setting the LED pin as output
}

//......................................Functions......................................

// Function for calculating the thermistor resistance through a voltage divider
float calculate_resistance(float v_th){               // v_th is the thermistor voltage
  float resistance = (v_in/v_th - 1)*R_2;             // Resistance calculation
  return resistance;                                  // Returning the resistance value         
}

// Function for calculating the thermistor temperature through the resistance value
float calculate_temperature(float r_th){              // r_th is the thermistor resistance
  float kelvin = 1 / (1/T_0 +  log(r_th/R_2)/beta);  // Kelvin temperature calculation
  float celsius = kelvin - KELVIN_ZERO;               // Celsius temperature calculation
  return celsius;                                     // Returning the temperature value         
}

// Function for changing the state of the LED according to the Thermistor temperature
void too_hot_to_handle(float t_th){                   // t_th is the thermistor temperature
  if (t_th > TARGET_TEMP)                             // If the temperature is higher than the target one 
    digitalWrite(LED_PIN, HIGH);                      // Turn the LED on
  else                                                // If the temperature is lower than the target temperature              
    digitalWrite(LED_PIN, LOW);                       // Turn the LED off
}

// Function for printing in the console 
void print_data(int rawInput, float voltage, float resistance, float temperature){ 
  Serial.print(rawInput);                             // Printing the raw input value
  Serial.print(';');                                  // Printing the separator 
  Serial.print(voltage);                              // Printing the voltage value
  Serial.print(';');
  Serial.print(resistance);                           // Printing the resistance value
  Serial.print(';');
  Serial.println(temperature);                        // Printing the temperature value
}

//.....................................Main Loop......................................
void loop() {
  // Reading the value from analog pin ANALOG_PIN
  int thermistorInput = analogRead(ANALOG_PIN);

  // Converting the analog input value to voltage
  float thermistorVoltage = thermistorInput * 5.0/1023;  

  // Calculating the thermistor resistance
  float thermistorResistance = calculate_resistance(thermistorVoltage);

  // Calculating the temperature value in Celsius 
  float thermistorTemperature = calculate_temperature(thermistorResistance);

  // Controlling the LED state according to the thermistor temperature
  too_hot_to_handle(thermistorTemperature);
  
  // Printing the data
  print_data(thermistorInput, thermistorVoltage, thermistorResistance, thermistorTemperature); 
  
  // Delaying the loop
  delay(DELAY);
}