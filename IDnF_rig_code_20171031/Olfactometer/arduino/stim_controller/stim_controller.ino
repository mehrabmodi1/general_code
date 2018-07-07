//Mehrab Modi, 20180707. This program receives 5 numbers via serial, denoted by '<' and '>'. It then uses these to trigger LED or electrical stimulation. The 5 numbers are
// select LED/elec (0/1), initial delay from trigger in ms, total stim dur in ms, freq in Hz and duty cycle as percent.

const byte numChars = 32;
char receivedChars[numChars];

boolean newData = false;
const int ledPin = 5;

void setup() {
    
    pinMode(ledPin, OUTPUT);
    Serial.begin(9600);
    Serial.println("<Arduino is ready>");
}

void loop() 
{   
    int param_n = 0;    //this variable encodes the number of newly received param values
    int LED_elec = 0;
    int init_delay_ms = 0;
    int duration_ms = 0;
    int freq_hz = 0;
    int duty_cycle_pc = 0;
    
    //reading in the 5 stimulus control parameters 
    while (param_n < 5)
    {
          recvWithStartEndMarkers();
          
          if (newData == true && param_n == 0) {LED_elec = atoi(receivedChars); param_n = param_n + 1; showNewData(); Serial.println(param_n);}
          if (newData == true && param_n == 1) {init_delay_ms = atoi(receivedChars); param_n = param_n + 1; showNewData(); Serial.println(param_n);}
          if (newData == true && param_n == 2) {duration_ms = atoi(receivedChars); param_n = param_n + 1; showNewData(); Serial.println(param_n);}
          if (newData == true && param_n == 3) {freq_hz = atoi(receivedChars); param_n = param_n + 1; showNewData(); Serial.println(param_n);}
          if (newData == true && param_n == 4) {duty_cycle_pc = atoi(receivedChars); param_n = param_n + 1; showNewData(); Serial.println(param_n);}      

    }

    //doing stuff with read parameters
    
    if (param_n == 5)
    {
          int cyc_dur = 1000/freq_hz;
          Serial.print("cyc_dur ");
          Serial.println(cyc_dur);

          float duty_cycle = (float)duty_cycle_pc/100;
          float on_dur = (float)cyc_dur*duty_cycle;
          Serial.print("on_dur ");
          Serial.println(on_dur);

          float off_dur = (float) cyc_dur - on_dur;
          Serial.print("off_dur ");
          Serial.println(off_dur);
          
          int n_pulses = duration_ms/cyc_dur;
          Serial.print("n_pulses");
          Serial.println(n_pulses);

          //delivering stimulus pulses
          int pulse_n = 0;

          //insert trigger wait loop here
          
          delay(init_delay_ms);
          while (pulse_n < n_pulses)
          {
                digitalWrite(ledPin, HIGH);
                delay(on_dur);
                digitalWrite(ledPin, LOW);
                delay(off_dur);
                pulse_n = pulse_n + 1;
          }
          //This resets the state to waiting for new data. 
          param_n = 1;
          
    }

}



//SERIAL READ FUNCTIONS

void recvWithStartEndMarkers() {
    static boolean recvInProgress = false;
    static byte ndx = 0;
    char startMarker = '<';
    char endMarker = '>';
    char rc;
 
    while (Serial.available() > 0 && newData == false) {
        rc = Serial.read();

        if (recvInProgress == true) {
            if (rc != endMarker) {
                receivedChars[ndx] = rc;
                ndx++;
                if (ndx >= numChars) {
                    ndx = numChars - 1;
                }
            }
            else {
                receivedChars[ndx] = '\0'; // terminate the string
                recvInProgress = false;
                ndx = 0;
                newData = true;
                
            }
        }

        else if (rc == startMarker) {
            recvInProgress = true;
        }
    }
}

void showNewData() {
    if (newData == true) {
        Serial.print("This just in ... ");
        Serial.println(receivedChars);
        newData = false;
        
    }
}


