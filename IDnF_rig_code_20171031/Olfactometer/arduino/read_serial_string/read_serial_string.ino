/*
  Reading a serial ASCII-encoded string.

 This sketch demonstrates the Serial parseInt() function.
 It looks for an ASCII string of comma-separated values.
 It parses them into ints, and uses those to fade an RGB LED.

 Circuit: Common-Cathode RGB LED wired like so:
 * Red anode: digital pin 3
 * Green anode: digital pin 5
 * Blue anode: digital pin 6
 * Cathode : GND

 created 13 Apr 2012
 by Tom Igoe
 
 modified 14 Mar 2016
 by Arturo Guadalupi

 This example code is in the public domain.
 */

// pins for the LEDs:
const int redPin = 2;

void setup() {
  // initialize serial:
  Serial.begin(9600);
  // make the pins outputs:
  pinMode(redPin, OUTPUT);
 
}

void loop() 
{
  // if there's any serial available, read it:
  while (Serial.available() > 0) 
  {

    // look for the next valid integer in the incoming serial stream:
    int init_delay = Serial.parseInt();
    // do it again:
    int duration = Serial.parseInt();
    // do it again:
    int freq = Serial.parseInt();
    // do it again:
    int duty_cycle = Serial.parseInt();
    
    // look for the newline. That's the end of your
    // sentence:
    if (Serial.read() == '\n') 
    {
      int n_pulses;
      int cyc_duration;
      int on_duration;
      int off_duration;
      n_pulses = duration*freq;

      // triggering LED or elec:
      digitalWrite(redPin, HIGH);
      delay(duration);
      digitalWrite(redPin, LOW);
      delay(freq);

      
      // print the three numbers in one string as hexadecimal:
      Serial.print(init_delay);
      Serial.print(duration);
      Serial.println(freq);
    }
  }
}








