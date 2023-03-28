#include <math.h>
#define limit 2.0
#define p_out 8




void setup() {
  Serial.begin(9600);
  pinMode(p_out, OUTPUT);
}

void loop() {
  int input = analogRead(A0);                             //get the info from pin A0(info)
  float v = 5.0 * input / 1023.0;                         //make the input voltage
  bool lamp = 0;                                        //boolean variable if led is on or off

  //If the voltage is better than the limit the led will blink and the lamp variavle will be true
  if(v>limit){
    digitalWrite(p_out,HIGH);
    lamp =1;
  }else{
    digitalWrite(p_out,LOW);
    lamp =0;
  }
  //Print the values

  Serial.print(input);         //Print the input
  Serial.print(';');
  Serial.print(v);            //Print the voltage
  Serial.print(';');
  Serial.print(lamp);            //Print the resistance
  Serial.print(';');
  Serial.print('\n');         //change line in order to seperate each measurement


  delay(500);         //delay of 0.5s
}