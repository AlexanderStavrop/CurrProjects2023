#include <math.h>             // Math library for log function
// Photoresistor constants
# define threshold 4.8       // Photoresistor threshold value
// I/O constants
#define ANALOG_PIN A0         // Thermistor analog pin
#define LED_PIN 8             // LED digital pin
// Experiment constants
#define DELAY 500             // Delay between measurements


//.......................................Setup........................................
void setup() {
  Serial.begin(9600);         // Starting the serial communication  
  pinMode(LED_PIN, OUTPUT);   // Setting the LED pin as output
}

//......................................Functions......................................

// Function for changing the state of the LED according to the Thermistor temperature
bool blackHole(float volt){                       // t_th is the thermistor temperature
  if (volt > threshold){                             // If the temperature is higher than the target one 
    digitalWrite(LED_PIN, HIGH);                      // Turn the LED on
    return 1;
  }
  else{                                                // If the temperature is lower than the target temperature              
    digitalWrite(LED_PIN, LOW);                       // Turn the LED off
    return 0;
  }
}

// Function for printing in the console 
void print_data(int rawInput, float voltage, bool state){ 
  Serial.print(rawInput);                             // Printing the raw input value
  Serial.print(';');                                  // Printing the separator 
  Serial.print(voltage);                              // Printing the voltage value
  Serial.print(';');
  Serial.println(state);                           // Printing the resistance value
}

//.....................................Main Loop......................................
void loop() {
  // Reading the value from analog pin ANALOG_PIN
  int photoresistorInput = analogRead(ANALOG_PIN);

  // Converting the analog input value to voltage
  float thermistorVoltage = photoresistorInput * 5.0/1023;  

  // Controlling the LED state according to the thermistor temperature
  bool LedState = blackHole(thermistorVoltage);
  
  // Printing the data
  print_data(photoresistorInput, thermistorVoltage, LedState); 
  
  // Delaying the loop
  delay(DELAY);
}