function [od_pulse_frames] = compute_pulse_frames_train(curr_train, frame_time, stim_latency)
%Syntax:[od_pulse_frames] = compute_pulse_frames_train(curr_train, frame_time, stim_latency)
%This function computes a sequence of pairs of frame numbers for on and
%off times of odor pulses from a train of odor on and odor off delays. 
%This can be used as an input to add_stim_bar, for example.

valve_odor_delay = load('C:\Data\Data\Analysed_data\valve_odor_delay.mat');
valve_odor_delay = valve_odor_delay.valve_odor_delay;    %in ms
valve_odor_delay = valve_odor_delay./1000; %in s

% %testing variables
% stim_latency = 25;
% frame_time = 99.93;     %in ms 
% curr_train =         [0,    8.2222; ...
%                11.4785,    6.6357; ...
%                 2.8376,    6.2957; ...
%                 3.6396,    1.8006; ...
%                 7.5791,    1.1169; ...
%                13.0986,   12.2855; ...
%                 6.6676,    1.8153; ...
%                 8.5403,    6.8911; ...
%                 2.7026,   18.3931];
            
stim_points = curr_train;
   
for pulse_n = 1:size(stim_points, 1)

    if pulse_n > 1
        stim_points(pulse_n, 1) = stim_points((pulse_n - 1), 2) + stim_points(pulse_n, 1);
        stim_points(pulse_n, 2) = stim_points(pulse_n, 1) + stim_points(pulse_n, 2);
    else
    end
end

od_pulse_frames = round(( (stim_points.*1) + (stim_latency.*1) + valve_odor_delay)./(frame_time));
