function [] = cell_data_quality_control(dff_data_mat, stim_mat, stim_mat_simple, sig_cell_mat, manual_inspec)

odor_list = unique(stim_mat_simple(:, 2));
n_trains = max(stim_mat_simple(:, 11));

global color_vec
global greymap
global frame_time

for odor_n = 1:length(odor_list)
    odor_ni = odor_list(odor_n);
    for train_n = 1:n_trains
        curr_od_trs = find(stim_mat_simple(:, 2) == odor_ni);
        curr_trn_trs = find(stim_mat_simple(:, 11) == train_n);
        curr_trs = intersect(curr_od_trs, curr_trn_trs);           %current trs
        curr_train = stim_mat(curr_trs(1)).rand_trains;
        od_pulse_frames = compute_pulse_frames_train(curr_train, frame_time, stim_mat(curr_trs(1)).stim_latency);
        ave_mat = mean(dff_data_mat(:, :, curr_trs), 3, 'omitnan');
        curr_sig_cells = sig_cell_mat(:, odor_ni);
        
        
        %making plots for manual inspection only - to check if sig cells have been found or left out
        if manual_inspec == 1
            curr_sig_cells1 = find(curr_sig_cells == 1);
            curr_sig_cells0 = find(curr_sig_cells == 0);
            figure(1)
            imagesc(ave_mat(:, curr_sig_cells1)', [0, 2])
            title('significant responders, all trials averaged')
            fig_wrapup(1)
            add_stim_bar(1, od_pulse_frames, color_vec(odor_ni, :));
            keyboard
            figure(2)
            imagesc(ave_mat(:, curr_sig_cells0)', [0, 2])
            title('not significant responders, all trials averaged')
            colormap(greymap)
        else
        end
        curr_sig_cells = find(curr_sig_cells == 1);    
        
        %checking if responses in first-last quarter of trials are similar to each other and to overall mean.
        q_n = floor(length(curr_trs)./4);
        
        if q_n > 0
            ave_mat_early = mean(dff_data_mat(:, curr_sig_cells, curr_trs(1:q_n)), 3, 'omitnan');
            ave_mat_late = mean(dff_data_mat(:, curr_sig_cells, curr_trs(((3*q_n) + 1): end)), 3, 'omitnan');

            curr_sig_cells = sig_cell_mat(:, odor_ni);
            curr_sig_cells = find(curr_sig_cells == 1);

            if manual_inspec == 1
                figure(3)
                keyboard
                imagesc(ave_mat_early(:, curr_sig_cells)', [0, 1.2])
                title('ave resp traces for first quarter of trials')
                colormap(greymap)
                figure(4)
                imagesc(ave_mat_late(:, curr_sig_cells)', [0, 1.2])
                title('ave resp traces for last quarter of trials')
                colormap(greymap)
            else
            end
            %computing difference between matrices
            keyboard
            
            
            
        else
            disp('too few reps to analyse by quarters of reps... skipping.')
        end
        keyboard

    end
end