clear all
close all

direc_list_path = 'C:\Data\Data\Raw_data\dataset_lists\dataset_list_Cafree_Ringers_expt.xls';
[del, direc_list] = xlsread(direc_list_path, 1);
n_direcs = size(direc_list, 1); 

for direc_n = 1:n_direcs
    curr_dir = direc_list{direc_n, 1};
    curr_dir = [curr_dir, '\'];
    remove_small_tifs(curr_dir);
    dir_contents = dir_date_sorted(curr_dir, '*.tif');
    [stim_mat, stim_mat_simple, column_heads] = load_params_trains(curr_dir, []);
    
    n_trials = size(stim_mat_simple, 1);
    
    for trial_n = 1:n_trials
        stack_obj = ScanImageTiffReader([curr_dir, dir_contents(trial_n).name]);
        stack = stack_obj.data();
        stack = permute(stack,[2 1 3]);
        n_frames = size(stack, 3);
        [frame_time, zoom, n_chans] = SI_tif_info(stack_obj);
        stim_fr = floor(stim_mat_simple(trial_n, 7)./frame_time);
        stim_end_fr = floor( (stim_mat_simple(trial_n, 7) + stim_mat_simple(trial_n, 3))./frame_time);
        
        highlight_resp_pix(1, stack, [stim_fr, stim_end_fr], 0.05, frame_time);
        
        
        
        
        
        keyboard
        
    end
    
    
end