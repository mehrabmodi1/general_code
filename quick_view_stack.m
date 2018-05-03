clear all
close all

curr_dir = 'C:\Data\Data\Raw_data\20180416\OK107_opGC6f\';

remove_small_tifs(curr_dir);
dir_contents = dir_date_sorted(curr_dir, '*.tif');
[stim_mat, stim_mat_simple, column_heads] = load_params_trains(curr_dir, []);
n_trials = size(stim_mat_simple, 1);

imagesc(stim_mat_simple, [0, 12])

trial_n = input('Enter trial number to view:');

stack_obj = ScanImageTiffReader([curr_dir, dir_contents(trial_n).name]);
stack = stack_obj.data();
stack = permute(stack,[2 1 3]);
n_frames = size(stack, 3);
[frame_time, zoom, n_chans] = SI_tif_info(stack_obj);

stim_fr = floor(stim_mat_simple(trial_n, 7)./frame_time);
stim_end_fr = floor( (stim_mat_simple(trial_n, 7) + stim_mat_simple(trial_n, 3))./frame_time);

highlight_resp_pix(1, stack, [stim_fr, stim_end_fr], 1, frame_time, 0);

fig_wrapup(1)