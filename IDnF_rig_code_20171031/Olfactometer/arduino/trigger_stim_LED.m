function [] = trigger_stim_LED(on_dur, off_dur, n_pulses) 

trig_ard_port = 'COM10';    %this is the serial port number for the trigger control arduino
trig_line = 'D3';           %this is the digital line number connected to Scanimage to trigger acqn.

iscon = check_ser_connected(trig_ard_port);

%connecting to arduino if necessary
persistent trig_ard
if iscon == 0
   trig_ard= arduino(trig_ard_port, 'Mega2560');
else
end    

try
    %loop for multiple pulses
    for pulse_n = 1:n_pulses
        
        writeDigitalPin(trig_ard, trig_line, 1);    %on part of the pulse
        pause(on_dur)                                      
        writeDigitalPin(trig_ard, trig_line, 0);    %off part of the pulse
        pause(off_dur)
        
    end
    
catch
    keyboard
end
