clear all
close all

direc = 'C:\Data\Data\Raw_data\20200820\fly4_vGlut_GC6fDrep\';
save_path = 'C:\Data\Data\Analysed_data\Analysis_results\vGlut_GCaMPDrep2\';

cd(direc);
remove_small_tifs(direc);
tif_list = dir('*.tif');

[stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig] = load_params_trains_modular(direc, []); 

od_olf1_col_n = find_stim_mat_simple_col('odor_n', column_heads);           %identifying relevant column number in stim_mat_simple
dur_olf1_col_n = find_stim_mat_simple_col('duration', column_heads);        %identifying relevant column number in stim_mat_simple
od_list = unique(stim_mat_simple(:, od_olf1_col_n));
od_durs = unique(stim_mat_simple(:, dur_olf1_col_n));
od_durs(isnan(od_durs)) = [];

if exist([save_path, 'ave_stacks.mat']) ~= 2
    ave_stacks = [];
    for od_ni = 1:length(od_list) 
        od_n = od_list(od_ni);
        curr_trs = find(stim_mat_simple(:, od_olf1_col_n) == od_n);
        for trial_ni = 1:length(curr_trs)
            trial_n = curr_trs(trial_ni);

            stack_obj = ScanImageTiffReader([direc, tif_list(trial_n).name]);
            [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);
            curr_stack = stack_obj.data();
            curr_stack = permute(curr_stack,[2 1 3]);

            if trial_ni == 1
                ave_stack = curr_stack;
            else
                ave_stack = pad_n_concatenate(ave_stack, curr_stack, 4, nan);
            end

        end    
        ave_stack = mean(ave_stack, 4, 'omitnan');
        ave_stacks = pad_n_concatenate(ave_stacks, ave_stack, 4, nan);
    end
    save([save_path, 'ave_stacks.mat'], 'ave_stacks');
    
else
    ave_stacks = load([save_path, 'ave_stacks.mat']);
    ave_stacks = ave_stacks.ave_stacks;
    stack_obj = ScanImageTiffReader([direc, tif_list(1).name]);
    [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);
end
keyboard
stim_frs = compute_stim_frs_modular(stim_mat, 1, frame_time);
stim_frs = stim_frs{1}; 
highlight_resp_pix(1, squeeze(ave_stacks(:, :, :, 1)), stim_frs, n_chans, frame_time, 0);


keyboard
