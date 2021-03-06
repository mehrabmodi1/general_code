function [stim_mat] = load_params_res(direc, n_trials_t)
%sytax: [stim_mat] = load_params_res(direc, n_trials_t)

% direc = 'D:\Data\CSHL\Analysed_data\Resonant\Results\20170509\odor_dur_trials\1\';
% n_trials_t = 38;                              %n_trials according to the loaded raw traces file

tif_datenums = load([direc, '\tif_datenums.mat']);
tif_datenums = tif_datenums.datenum_vec;

%identifying newest params file and loading it
prev_direc = pwd;
cd(direc)
dir_contents = dir('params*.mat' );
n_files = size(dir_contents, 1);

date_nums = zeros(n_files, 1) + nan;
for par_file_n = 1:n_files
    date_nums(par_file_n) = dir_contents(par_file_n).datenum;
end

[del, last_filen] = max(date_nums);
par_filename = dir_contents(last_filen).name;
params = load([direc, par_filename]);
params = params.params;
n_trials_p = size(params.timestamp, 2);        %n_trials according to the param file

%n_trials_p = 3;
n_trials = min([n_trials_t, n_trials_p]);
match_mat = zeros(n_trials_p, n_trials_t) + nan;

for trial_n_p = 1:n_trials_p
    curr_datenum_p = params.timestamp(trial_n_p);
    
    for trial_n_t = 1:n_trials_t
        curr_datenum_t = tif_datenums(trial_n_t);
        match_mat(trial_n_p, trial_n_t) = etime(datevec(curr_datenum_t), datevec(curr_datenum_p));      %calculating time elapsed from param time stamp to tif time stamp             
    end   
    
end

match_valsi = find(match_mat < 2);         %finding pairs of tif and param datenums less than 2 s apart
match_valsii = find(match_mat > -2);
match_valsi = intersect(match_valsi, match_valsii);
[par_num, tif_num] = ind2sub([n_trials_p, n_trials_t], match_valsi);        

% figure(1)
% imagesc(match_mat)
% xlabel('tif trial n');
% ylabel('param trial n');
% 
% figure(2)
% plot(tif_num, par_num, '.');
% xlabel('tif trial n')
% ylabel('param trial n');

n_matched_trials = length(par_num);

stim_mat_dat = zeros(n_matched_trials, 9);
column_heads = '[matched_tif_n, odor_n, duration, isi, n_odor_pulses, inter_pulse_interval, stim_latency, first_dilution, second_dilution]';
%loop to read param values from params file into stim_mat
for trial_n = 1:n_matched_trials
    curr_tr_p = par_num(trial_n);
    curr_tr_t = tif_num(trial_n);
    odor_n = params.odours(curr_tr_p);
    duration = params.duration(curr_tr_p);
    isi = params.isi(curr_tr_p);
    n_odor_pulses = params.n_od_pulses(curr_tr_p);
    inter_pulse_interval = params.inter_pulse_interval(curr_tr_p);
    stim_latency = params.stimLatency(curr_tr_p);
    first_dilution = params.firstDilution(curr_tr_p);
    second_dilution = params.secondDilution(curr_tr_p);
    stim_mat_dat(trial_n, :) = [curr_tr_t, odor_n, duration, isi, n_odor_pulses, inter_pulse_interval, stim_latency, first_dilution, second_dilution];
end

save([direc, '\stim_mat.mat'], 'stim_mat')

end