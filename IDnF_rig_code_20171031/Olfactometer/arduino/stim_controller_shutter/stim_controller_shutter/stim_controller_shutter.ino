//Mehrab Modi, 20180707. This program receives 5 numbers via serial, denoted by '<' and '>'. It then uses these to trigger LED or electrical stimulation. The 5 numbers are
// select LED/elec (0/1), initial delay from trigger in ms, total stim dur in ms, freq in Hz and duty cycle as percent. It also generates a warning signal whenever the stim 
//LED is about to come on so as to close the PMT shutter and protect it.

const byte numChars = 32;
char receivedChars[numChars];

boolean newData = false;
const int led_pin = 5;
const int elec_pin = 6;
const int trig_pin = 4;
const int led_warning_pin = 7;
const float led_warning_leadt = 20;

int trig_state = 0;
int param_n = 0;    //this variable encodes the number of newly received param values
int LED_elec = 0;
float init_delay_ms = 0;
float duration_ms = 0;
float freq_hz = 0;
float duty_cycle_pc = 0;
float on_dur_us = 0;

void setup() {
    
    pinMode(led_pin, OUTPUT);
    pinMode(elec_pin, OUTPUT);
    pinMode(trig_pin, INPUT);
    Serial.begin(9600);
    Serial.println("<Arduino is ready>");
}

void loop() 
{   
    
    
    //reading in the five stimulus control parameters 
    while (param_n < 5)
    {
          recvWithStartEndMarkers();
          
          if (newData == true && param_n == 0) {LED_elec = atoi(receivedChars); param_n = param_n + 1; showNewData(); Serial.println(param_n);}
          if (newData == true && param_n == 1) {init_delay_ms = atof(receivedChars); param_n = param_n + 1; showNewData(); Serial.println(param_n);}
          if (newData == true && param_n == 2) {duration_ms = atof(receivedChars); param_n = param_n + 1; showNewData(); Serial.println(param_n);}
          if (newData == true && param_n == 3) {freq_hz = atof(receivedChars); param_n = param_n + 1; showNewData(); Serial.println(param_n);}
          if (newData == true && param_n == 4) {duty_cycle_pc = atof(receivedChars); param_n = param_n + 1; showNewData(); Serial.println(param_n);}      

    }

    //beginning stimulus delivery only if all five parameters have been received
    if (param_n == 5)
    {     
          //computing stimulus delivery variables
          float cyc_dur = (float) 1000/(float)freq_hz;
          Serial.print("cyc_dur ");
          Serial.println(cyc_dur);

          float duty_cycle = (float)duty_cycle_pc/100;
          float on_dur = (float)cyc_dur*duty_cycle;
          Serial.print("on_dur ");
          Serial.println(on_dur);

          float off_dur = (float) cyc_dur - on_dur;
          Serial.print("off_dur ");
          Serial.println(off_dur);
          float off_dur_w = (float) off_dur - (float)led_warning_leadt;
          
          float n_pulses = (float)duration_ms/(float)cyc_dur;
          Serial.print("n_pulses");
          Serial.println(n_pulses);


          //delivering stimulus pulses
          
          //waiting for scan trigger
          Serial.print("waiting for trig..");
          while (trig_state == LOW) {trig_state = digitalRead(trig_pin);}

          //waiting for initial delay, but generating led warning if needed. This assumes that initial delay > led warning lead time
          float init_delay_w = (float) init_delay_ms - led_warning_leadt;
          delay(init_delay_w); //waiting for initial delay after scan trigger
          

          //delivering stim pulses       
          float pulse_n = 0;

          while (pulse_n < n_pulses)
          {
                if (LED_elec == 0) {digitalWrite(led_warning_pin, HIGH);}   //Turning on LED warning signal
                delay(led_warning_leadt);

                //Turning on LED or elec stim
                if (LED_elec == 0) {digitalWrite(led_pin, HIGH);}           
                if (LED_elec == 1) {digitalWrite(elec_pin, HIGH);}
                
                if (on_dur > 2)
                  {
                    delay(on_dur);    //pulse on duration
                  }
                else if (on_dur <= 2)
                  {
                    on_dur_us = on_dur * 1000;
                    delayMicroseconds(on_dur_us);
                  }
                
                //Turning off LED or elec stim
                if (LED_elec == 0) {digitalWrite(led_pin, LOW);}
                if (LED_elec == 1) {digitalWrite(elec_pin, LOW);}
                
                if (off_dur > led_warning_leadt) 
                        {digitalWrite(led_warning_pin, LOW);  
                          delay(off_dur_w);}                  //pulse off duration  //turning off led warning pin during inter pulse interval
                
                if (off_dur <= led_warning_leadt) {delay(off_dur);}         //if off dur is too short for LED warning, led_warning_pin is left high and off_dur is used.

                
                pulse_n = pulse_n + 1;
          }
          digitalWrite(led_warning_pin, LOW);
          
          //This resets the state to waiting for new parameter data and a new scan  trigger.
          param_n = 0;
          trig_state = 0;
          Serial.print("done with stim.");
         
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


