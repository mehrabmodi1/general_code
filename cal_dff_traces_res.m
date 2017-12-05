function [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, direc)
%This function takes the extracted raw traces in a 3-D matrix (frames,
%cells, trials) and the corress 2P data object and gives an output of dF/F
%traces in a matrix identical in size to raw_data_mat, with an extra dimension for odors.  

% direc = 'C:\Data\CSHL\20150227\expt\';
% 
% dataset = load([direc 'expt.mat']);
% dataset = dataset.data;
% raw_data_mat = load([direc 'expt_raw_traces.mat']);
% raw_data_mat = raw_data_mat.raw_data_mat;

PID_latency = 200;                                          %latency from odor onset to arrival at the end of the odor tube measured with a PID (in ms)

bad_tr_list = load([direc, '\bad_trial_list.mat']);
good_tr_list = bad_tr_list.bad_tr_list;                     %this is actually the list of good trials

n_frames = size(raw_data_mat, 1);                           %this is the maximum number of frames for any trial in this dataset
n_cells = size(raw_data_mat, 2);
n_trials = size(raw_data_mat, 3);

dff_data_mat = zeros(n_frames, n_cells, n_trials) + nan; %dff_data_mat can accommodate all possible odor numbers
dff_data_mat_f = zeros(n_frames, n_cells, n_trials) + nan; %dff_data_mat can accommodate all possible odor numbers


for trial_n = 1:n_trials
    stim_time = stim_mat(trial_n).stim_latency;
    
    stim_frame = floor(stim_time./frame_time) + PID_latency;
    odor_duration = stim_mat(trial_n).odor_duration;
    
    %calculating dF/F for all cells, for current trial
    raw_mat = squeeze(raw_data_mat(:, :, trial_n) );
    
    baseline_vec = mean(raw_mat(1:(stim_frame - 2), :), 1, 'omitnan');          %vector of F baselines for all cells
    baseline_mat = repmat(baseline_vec, n_frames, 1);
    dff_mat = (raw_mat - baseline_mat)./baseline_mat;
    b_vec = zeros(1, round(filt_time./frame_time)) + 1;
    a_vec = round(1000./frame_time);
    dff_mat_f = filter(b_vec, a_vec, dff_mat, [], 1);
        
    %checking if this trial was identified as bad and skipping if so.
    if isempty(intersect(good_tr_list, trial_n)) == 1
        continue            %leaves this trial as NaNs in dff_data_mat
    else
        try
            dff_data_mat(:, :, trial_n) = dff_mat;
            dff_data_mat_f(:, :, trial_n) = dff_mat_f;
        catch
        	keyboard
        end
    end
        
    
end

end