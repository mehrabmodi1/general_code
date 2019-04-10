/*
  modified from Button Example
 */

// constants won't change. They're used here to
// set pin numbers:
const int ScanImagePMTPinG = 2;     // Scanimage control signal in (green)  
const int ScanImagePMTPinR = 3;     // Scanimage control signal in (red)
const int LEDWarningPin = 4;        // LED warning signal in from stim LED controller
const int shutterOutPinG =  6;      // signal out to shutter controller (green)
const int shutterOutPinR =  7;      // signal out to shutter controller (red)



void setup() {
  // initialize the shutter control pins as an outputs:
  pinMode(shutterOutPinG, OUTPUT);
  pinMode(shutterOutPinR, OUTPUT);
  // initialize the SI PMT inputs and stim LED warning pins as an inputs:
  pinMode(ScanImagePMTPinG, INPUT);
  pinMode(ScanImagePMTPinR, INPUT);
  pinMode(LEDWarningPin, INPUT);
}

void loop() {
  // read the state of the pushbutton value:
 bool SIstateG = digitalRead(ScanImagePMTPinG);
 bool SIstateR = digitalRead(ScanImagePMTPinR);
 bool LEDstate = digitalRead(LEDWarningPin);
  
  // opening G_shutter only if LEDstate is low and G_SIstate is high
  if (SIstateG == HIGH && LEDstate == LOW) {
    // opening shutter:
    digitalWrite(shutterOutPinG, HIGH);
  } else {
    digitalWrite(shutterOutPinG, LOW);
  }

  // opening R_shutter only if LEDstate is low and G_SIstate is high
  if (SIstateR == HIGH && LEDstate == LOW) {
    // opening shutter:
    digitalWrite(shutterOutPinR, HIGH);
  } else {
    digitalWrite(shutterOutPinR, LOW);
  }
}
