function [] = trigger_scan(write_val, write_val_olf2) 

trig_ard_port = 'COM10';    %this is the serial port number for the trigger control arduino
trig_line = 'D2';           %this is the digital line number connected to Scanimage to trigger acqn.
trig_line_olf2 = 'D4';

iscon = check_ser_connected(trig_ard_port);
if nargin == 1
    write_val_olf2 = write_val;
else
end

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
        
        if strcmp(write_val_olf2, 'high') == 1
            writeDigitalPin(trig_ard, trig_line_olf2, 1);
        elseif strcmp(write_val_olf2, 'low') == 1
            writeDigitalPin(trig_ard, trig_line_olf2, 0);
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
        
        if write_val_olf2 == 1
            writeDigitalPin(trig_ard, trig_line_olf2, 1);

        elseif write_val_olf2 == 0
            writeDigitalPin(trig_ard, trig_line_olf2, 0);
        else
            error('acceptable inputs are 0, 1, ''low'' and ''high''. ')
        end


    end
catch
    keyboard
end