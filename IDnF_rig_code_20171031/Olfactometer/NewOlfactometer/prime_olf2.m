function [] = prime_olf2(od_list)

if isempty(od_list) == 1
    od_list = 1:4;
else
end

n_pulses = length(od_list).*2;

pulse_train = zeros(n_pulses, 2) + 5;
pulse_train(:, 2) = 8;
od_vec = repmat(od_list, 1, 2);
od_vec = reshape(od_vec, 1, []);

olf_arduino_serial_comm(0, pulse_train, od_vec, 1);
pause(0.5)

trigger_scan(1);
pause(1)
trigger_scan(0);
pause(size(pulse_train, 1).*10);

sleep_olf2;