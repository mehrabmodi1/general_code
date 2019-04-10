function [dff_data_mat, stim_mat] = cal_dff_traces(raw_data_mat, dataset);
%This function takes the extracted raw traces in a 3-D matrix (frames,
%cells, trials) and the corress 2P data object and gives an output of dF/F
%traces in a matrix identical in size to raw_data_mat, with an extra dimension for odors.  

% direc = 'C:\Data\CSHL\20150227\expt\';
% 
% dataset = load([direc 'expt.mat']);
% dataset = dataset.data;
% raw_data_mat = load([direc 'expt_raw_traces.mat']);
% raw_data_mat = raw_data_mat.raw_data_mat;


stim_time = dataset(1).stim.stimLatency.*1000;              %stimulus onset time in ms
stim_time = stim_time + 625;                                %added delay from valve opening to odor at pipe outlet
frame_time = dataset(1).info.framePeriod .* 1000;           %frame time in ms
n_frames = size(raw_data_mat, 1);
n_cells = size(raw_data_mat, 2);
n_trials = size(raw_data_mat, 3);
stim_frame = floor(stim_time./frame_time);                  %frame no at which odor reached fly

clear frame_time


dff_data_mat = zeros(n_frames, n_cells, n_trials, 8) + nan; %dff_data_mat can accommodate all possible odor numbers
stim_mat = zeros(n_trials, 4);

try
    treatment_switch_trial = dataset(1).stim.switch_valve_trial;
    if isempty(treatment_switch_trial) == 1
        treatment_switch_trial = 1;
    else
    end
catch
    treatment_switch_trial = 1;
end

for trial_n = 1:n_trials
    odor_n = dataset(trial_n).stim.odours;
    
    %checking if elec stim was delivered on this trial
    elec_odors = dataset(trial_n).stim.elec_odours;
    if isempty(elec_odors) == 1
        elec_stim = 0;
    else
        elec_odors = str2num(elec_odors);
        if isempty(intersect(elec_odors, odor_n)) == 1
            elec_stim = 0;
        elseif isempty(intersect(elec_odors, odor_n)) == 0
            elec_stim = 1;
        else
        end
    end
    
    %checking if led stim was delivered on this trial
    led_odors = dataset(trial_n).stim.led_odours;
    if isempty(led_odors) == 1
        led_stim = 0;
    else
        led_odours = str2num(led_odors);
        if isempty(intersect(led_odors, odor_n)) == 1
            led_stim = 0;
        elseif isempty(intersect(led_odors, odor_n)) == 0
            led_stim = 1;
        else
        end
    end
    
    
    if trial_n > treatment_switch_trial
        treatment_state = 1;
    else 
        treatment_state = 0;
    end
    stim_vec = [odor_n, elec_stim, led_stim, treatment_state];
    stim_mat(trial_n, :) = stim_vec;                    %saving details of stimuli delivered
    
    %calculating dF/F for all cells, for current trial
    raw_mat = squeeze(raw_data_mat(:, :, trial_n) );
    
    baseline_vec = nanmean(raw_mat(1:(stim_frame - 2), :), 1);
    baseline_mat = repmat(baseline_vec, n_frames, 1);
    dff_mat = (raw_mat - baseline_mat)./baseline_mat;
    
    dff_data_mat(:, :, trial_n, odor_n) = dff_mat;

end

clear stim_vec
clear elec_odors
clear led_odors
clear dataset
clear n_cells
clear n_frames
clear n_trials
clear odor_n
clear raw_data_mat
clear stim_frame
clear stim_time
clear trial_n
clear baseline_mat
clear baseline_vec
clear dff_mat
clear elec_stim
clear led_stim
clear raw_mat


end