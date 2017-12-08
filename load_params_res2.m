function [stim_mat_simple, column_heads, rand_trains] = load_params_res2(params)
%This function loads in a param file as per 12/8/2017 and parses it into a
%stim_mat structure and a stim_mat_simple matrix whose column headers are
%also specified as strings.

%test vars
% clear all
% close all
% direc = 'D:\Data\Janelia\resonant\20171206\olf_testing2\';
% params = load([direc, 'params.mat']);
% params = params.params_mat;

n_trials = size(params, 2);
column_heads = [{'odor_n'}, {'duration'}, {'isi'}, {'n_odor_pulses'}, {'inter_pulse_interval'}, {'stim_latency'}, {'first_dilution'}, {'second_dilution'},...
                    {'post_od_scan_dur'}, {'rand_train_n'}];

stim_mat_simple = zeros(n_trials, 10) + nan;
rand_trains = [];
for trial_n = 1:n_trials
    param_vec = [params(trial_n).odours, params(trial_n).duration, params(trial_n).isi,...
        params(trial_n).n_od_pulses, params(trial_n).inter_pulse_interval, params(trial_n).stimLatency,...
        params(trial_n).firstDilution, params(trial_n).secondDilution, params(trial_n).post_od_scan_dur,...
        params(trial_n).rand_train_n];
    
    stim_mat_simple(trial_n, :) = param_vec;
    
    
    rand_train_n = params(trial_n).rand_train_n;
    rand_trains{rand_train_n, 1} = params(trial_n).rand_train;    
end




%end