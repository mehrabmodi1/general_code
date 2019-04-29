function stim_arduino_serial_comm(LED_elec, init_delay_ms, duration_ms, freq_hz, duty_cyc_percent)

stim_ard = serial('COM13');
set(stim_ard,'BaudRate',9600);
fopen(stim_ard);

pause(2);       %waiting for matlab to open the serial comm port

%communicating stimulus parameters to stim arduino
fprintf(stim_ard, '%s', ['<', num2str(LED_elec), '>']);              %specifying LED or elec stim (0 or 1 resp.)
fprintf(stim_ard, '%s', ['<', num2str(init_delay_ms), '>']);         %initial delay after scan trigger onset in ms
fprintf(stim_ard, '%s', ['<', num2str(duration_ms), '>']);           %total stim duration in ms 
fprintf(stim_ard, '%s', ['<', num2str(freq_hz), '>']);               %frequency of stim pulses in Hz
fprintf(stim_ard, '%s', ['<', num2str(duty_cyc_percent), '>']);      %duty cycle of on part of stim pulse in percent

%closing, clearing stim arduino comm port
fclose(stim_ard);
delete(stim_ard);
clear stim_ard