function olf_arduino_serial_comm(mid_trial, pulse_train, odor_vec, initial_delay)
port_number = 13;

%making sure port is closed
close_serial_port(port_number);
port_name = ['COM', num2str(port_number)];
olf2_ard = serial(port_name);
set(olf2_ard,'BaudRate',9600);
fopen(olf2_ard);
pause(2);       %waiting for matlab to open the serial comm port

%figuring out if this is the first or last trial
if mid_trial == 0       %case when this is the first trial of a dataset
    final_trial = 0;
    fprintf(olf2_ard, '%s', ['<', 1, '>']);          %waking up olfactometer
    pause(0.5);
elseif mid_trial == 2   %case when this is the last trial of a dataset
    final_trial = 1;
elseif mid_trial == 1   %case when this is a middle trial (nothing special)
    final_trial = 0;
end

%converting s into ms
initial_delay = initial_delay.*1000;
pulse_train = pulse_train.*1000;
n_pulses = size(pulse_train, 1);

%communicating stimulus parameters to stim arduino
try
    fprintf(olf2_ard, '%s', ['<', num2str(final_trial), '>']);           %specifying whether the olfactometer should go to sleep after this trial
catch
    keyboard
end
pause(0.05);
fprintf(olf2_ard, '%s', ['<', num2str(n_pulses), '>']);              %initial delay after scan trigger onset in ms
pause(0.05);
%loop to communicate each pulse's parameters

for pulse_n = 1:n_pulses
    fprintf(olf2_ard, '%s', ['<', num2str(pulse_train(pulse_n, 2) ), '>']);     %odor on time in ms 
    pause(0.05);
end
for pulse_n = 1:n_pulses
    fprintf(olf2_ard, '%s', ['<', num2str(pulse_train(pulse_n, 1)), '>']);      %odor_off time in ms
    pause(0.05);
end
for pulse_n = 1:n_pulses    
    fprintf(olf2_ard, '%s', ['<', num2str(odor_vec(pulse_n)), '>']);            %odor_n for current pulse
    pause(0.05);
end

%communicating last parameter
fprintf(olf2_ard, '%s', ['<', num2str(initial_delay), '>']);             %LED control voltage to be delivered

%closing, clearing stim arduino comm port if this is the last trial
if mid_trial == 2
    fclose(olf2_ard);
    delete(olf2_ard);
    clear olf2_ard
else
end