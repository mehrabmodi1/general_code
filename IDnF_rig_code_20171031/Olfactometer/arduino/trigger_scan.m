function [] = trigger_scan(write_val) 

trig_ard_port = 'COM10';    %this is the serial port number for the trigger control arduino
trig_line = 'D2';           %this is the digital line number connected to Scanimage to trigger acqn.

iscon = check_ser_connected(trig_ard_port);

%connecting to arduino if necessary
persistent trig_ard
if iscon == 0
   trig_ard= arduino(trig_ard_port, 'Mega2560');
else
end    

try
    if ischar(write_val) == 1

        if strcmp(write_val, 'high') == 1
            writeDigitalPin(trig_ard, trig_line, 1);
        elseif strcmp(write_val, 'low') == 1
            writeDigitalPin(trig_ard, trig_line, 0);
        else
        end

    else
        if write_val == 1
            writeDigitalPin(trig_ard, trig_line, 1);

        elseif write_val == 0
            writeDigitalPin(trig_ard, trig_line, 0);
        else
            error('acceptable inputs are 0, 1, ''low'' and ''high''. ')
        end

    end
catch
    keyboard
end