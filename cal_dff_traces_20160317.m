function [dff_data_mat, stim_mat, prot_switch_trials] = cal_dff_traces_20160317(raw_data_mat, dataset, list_direc)
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
stim_time = stim_time + 450;                                %added delay from valve opening to odor at pipe outlet
frame_time = dataset(1).info.framePeriod .* 1000;           %frame time in ms
n_frames = size(raw_data_mat, 1);
n_cells = size(raw_data_mat, 2);
n_trials = size(raw_data_mat, 3);
stim_frame = floor(stim_time./frame_time);                  %frame no at which odor reached fly

clear frame_time


dff_data_mat = zeros(n_frames, n_cells, n_trials, 8) + nan; %dff_data_mat can accommodate all possible odor numbers
stim_mat = zeros(n_trials, 4) + nan;


%extracting protocol trial blocks in dataset and LED, elec odours for these
%blocks
prot_switch_trials = dataset(1).stim.prot_switch_trials;
del = findstr(list_direc, 'Toshi_KC_DN_Led_20150708');
if isempty(del) == 0
    prot_switch_trials = [10; 12; 36];
else
    prot_switch_trials = str2num(prot_switch_trials);
end



led_odors_mat = dataset(1).stim.led_odours;
elec_odors_mat = dataset(1).stim.elec_odours;
curr_block_n = 1;

if isempty(led_odors_mat) == 0
    led_odors = str2num(led_odors_mat(1, :));
else
    led_odors = [];
end

if isempty(elec_odors_mat) == 0
    elec_odors = str2num(elec_odors_mat(1, :));
else
    elec_odors = [];
end

for trial_n = 1:n_trials
    
    odor_n = dataset(trial_n).stim.odours;
    if isnan(odor_n) == 1
        continue
    else
    end
    odor_duration = dataset(trial_n).stim.duration;
        
    %checking if elec stim was delivered on this trial
    if isempty(elec_odors) == 1
        elec_stim = 0;
    else
        
        if isempty(intersect(elec_odors, odor_n)) == 1
            elec_stim = 0;
        elseif isempty(intersect(elec_odors, odor_n)) == 0
            elec_stim = 1;
        else
        end
    end
    
    %checking if led stim was delivered on this trial
    if isempty(led_odors) == 1
        led_stim = 0;
    else
        
        if isempty(intersect(led_odors, odor_n)) == 1
            led_stim = 0;
        elseif isempty(intersect(led_odors, odor_n)) == 0
            led_stim = 1;
            
        else
        end
    end
    
    
    stim_vec = [odor_n, odor_duration, elec_stim, led_stim];

    stim_mat(trial_n, :) = stim_vec;                    %saving details of stimuli delivered
    
    
    %calculating dF/F for all cells, for current trial
    raw_mat = squeeze(raw_data_mat(:, :, trial_n) );
    
    baseline_vec = nanmean(raw_mat(1:(stim_frame - 2), :), 1);          %vector of F baselines for all cells
    baseline_mat = repmat(baseline_vec, n_frames, 1);
    dff_mat = (raw_mat - baseline_mat)./baseline_mat;
    
    %checking if this trial was identified as bad and skipping if so.
    if dataset(trial_n).process.good_trial == 1
        dff_data_mat(:, :, trial_n, odor_n) = nan;
        
    else
        try
            dff_data_mat(:, :, trial_n, odor_n) = dff_mat;
        catch
        	keyboard
        end
    end
    
    %checking if end of protocol block has been reached and switching to
    %next protocol block if so.
    if curr_block_n < length(prot_switch_trials)
        if trial_n == prot_switch_trials(curr_block_n)
            curr_block_n = curr_block_n + 1;
            
            if isempty(led_odors_mat) == 0
                led_odors = str2num(led_odors_mat(curr_block_n, :));
            else
                led_odors = [];
            end
            del = find(led_odors > 8);
            led_odors(del) = [];
            
            if isempty(elec_odors_mat) == 0
                elec_odors = str2num(elec_odors_mat(curr_block_n, :));
            else
                elec_odors = [];
            end
            del = find(elec_odors > 8);
            elec_odors(del) = [];
            
        else
        end
    else
    end
    
    
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