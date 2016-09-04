function [response_mat, sd_mat, ave_resp_mat, ave_sd_mat] = find_pk_resp(dff_data_mat, stim_mat, frame_time, n_reps, stim_time)

win_width = 5;          %in n frames, the width of the window to calculate area under the dF/F curve.

n_cells = size(dff_data_mat, 2);
n_trials = size(dff_data_mat, 3);

odor_list = unique(stim_mat(:, 1));
n_odors = length(odor_list);
odor_dur_list = unique(stim_mat(:, 2));
n_odor_durs = length(odor_dur_list);

stim_frame = floor(stim_time./frame_time);
response_mat = zeros(n_cells, n_odors, n_odor_durs, n_reps) + nan;      %initialising response matrix
sd_mat = zeros(n_cells, n_odors, n_odor_durs, n_reps) + nan;
ave_resp_mat = zeros(n_cells, n_odors, n_odor_durs) + nan;      %initialising response matrix
ave_sd_mat = zeros(n_cells, n_odors, n_odor_durs) + nan;

for odor_dur_n = 1:n_odor_durs
    curr_od_dur = odor_dur_list(odor_dur_n);
    stim_end_fr = round((stim_time + curr_od_dur)./frame_time);
    curr_dur_trs = find(stim_mat(:, 2) == curr_od_dur);
        
    for odor_n = 1:n_odors
        odor_ni = odor_list(odor_n);
        curr_od_trs = find(stim_mat(:, 1) == odor_ni);
        curr_trs = intersect(curr_od_trs, curr_dur_trs);       %list of reps (max of 5) for curr od at curr duration
        
        if isempty(curr_trs) == 1
            continue
        else
        end
        
        curr_traces = squeeze(dff_data_mat(:, :, curr_trs, odor_ni));
        curr_ave_traces = nanmean(curr_traces, 3);
        
        curr_traces_an = curr_traces(stim_frame:(stim_end_fr + 5), :, :);
        curr_ave_traces_an = curr_ave_traces(stim_frame:(stim_end_fr + 5), :);
        
        baseline_traces = curr_traces( (stim_frame - 21):(stim_frame - 1), :, :);
        base_sds = squeeze(std(baseline_traces, [], 1));
        
        baseline_ave_traces = curr_ave_traces( (stim_frame - 21):(stim_frame - 1), :, :);
        base_ave_sds = squeeze(std(baseline_ave_traces, [], 1));
        
        
        %computing 5-frame moving window averaged trace and finding peak
        curr_traces_f = tsmovavg_m(curr_traces_an, 's', 1);
        curr_ave_traces_f = tsmovavg_m(curr_ave_traces_an, 's', 1);
                       
        max_resps = max(curr_traces_f, [], 1);
        max_ave_resps = max(curr_ave_traces_f, [], 1);
        
        response_mat(:, odor_n, odor_dur_n, 1:length(curr_trs)) = max_resps;
        
        sd_mat(:, odor_n, odor_dur_n, 1:length(curr_trs)) = base_sds;
        
        ave_resp_mat(:, odor_n, odor_dur_n) = max_ave_resps;
        ave_sd_mat(:, odor_n, odor_dur_n) = base_ave_sds;
    end
end

end