function [] = cell_data_quality_control(dff_data_mat, stim_mat_simple, sig_cell_mat)

odor_list = unique(stim_mat_simple(:, 2));
n_trains = max(stim_mat_simple(:, 11));

for odor_n = 1:length(odor_list)
    odor_ni = odor_list(odor_n);
    for train_n = 1:n_trains
        curr_od_trs = find(stim_mat_simple(:, 2) == odor_ni);
        curr_trn_trs = find(stim_mat_simple(:, 11) == train_n);
        curr_trs = intersect(curr_od_trs, curr_trn_trs);           %current trs
        half_n = floor(length(curr_trs)./2);
        q_n = floor(length(curr_trs)./4);

        ave_mat = mean(dff_data_mat(:, :, curr_trs), 3, 'omitnan');
        curr_sig_cells1 = find(curr_sig_cells == 1);
        curr_sig_cells0 = find(curr_sig_cells == 0);
        figure(1)
        imagesc(ave_mat(:, curr_sig_cells1)', [0, 1.2])
        title('significant responders, all trials averaged')
        figure(2)
        imagesc(ave_mat(:, curr_sig_cells0)', [0, 1.2])
        title('not significant responders, all trials averaged')
        
        ave_mat_early = mean(dff_data_mat(:, :, curr_trs(1:q_n)), 3, 'omitnan');
        ave_mat_late = mean(dff_data_mat(:, :, curr_trs(((3*q_n) + 1): end)), 3, 'omitnan');

        curr_sig_cells = sig_cell_mat(:, odor_ni);
        curr_sig_cells = find(curr_sig_cells == 1);

        figure(3)
        imagesc(ave_mat_early(:, curr_sig_cells)', [0, 1.2])
        figure(4)
        imagesc(ave_mat_late(:, curr_sig_cells)', [0, 1.2])

        keyboard

    end
end