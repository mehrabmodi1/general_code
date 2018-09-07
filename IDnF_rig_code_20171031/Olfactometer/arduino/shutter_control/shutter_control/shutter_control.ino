/*
  modified from Button Example
 */

// constants won't change. They're used here to
// set pin numbers:
const int ScanImagePMTPin = 2;     // the number of the pushbutton pin
const int LEDWarningPin = 4;     // the number of the pushbutton pin
const int shutterOutPin =  6;      // the number of the LED pin

// variables will change:
int buttonState = 0;         // variable for reading the pushbutton status

void setup() {
  // initialize the shutter control pin as an output:
  pinMode(shutterOutPin, OUTPUT);
  // initialize the SI PMT input and stim LED warning pins as an inputs:
  pinMode(ScanImagePMTPin, INPUT);
  pinMode(LEDWarningPin, INPUT);
}

void loop() {
  // read the state of the pushbutton value:
 bool SIstate = digitalRead(ScanImagePMTPin);
 bool LEDstate = digitalRead(LEDWarningPin);
  
  // opening shutter only if LEDstate is low and SIstate is high
  if (SIstate == HIGH && LEDstate == LOW) {
    // opening shutter:
    digitalWrite(shutterOutPin, HIGH);
  } else {
    // turn LED off:
    digitalWrite(shutterOutPin, LOW);
  }
}
