function [] = trigger_stim_LED_simple(write_val) 

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
    
     writeDigitalPin(trig_ard, trig_line, write_val);    %on part of the pulse   %on part of the pulse
        
catch
    keyboard
end
