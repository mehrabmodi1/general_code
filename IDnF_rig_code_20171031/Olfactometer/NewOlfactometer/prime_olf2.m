function [] = prime_olf2(od_list)

if isempty(od_list) == 1
    od_list = 1:4;
else
end

n_pulses = length(od_list).*2;

pulse_train = zeros(n_pulses, 2) + 5;
od_vec = repmat(od_list, 1, 2);
od_vec = reshape(od_vec, 1, []);

olf_arduino_serial_comm(0, odor_train, odor_vec, 1);
pause(0.5)

trigger_scan(1);
pause(1)
trigger_scan(0);

sleep_olf2;