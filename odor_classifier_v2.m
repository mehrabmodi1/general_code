function [classifier_output] = odor_classifier(sim_data_mat, duration, integration_window)

n_cells = size(sim_data_mat, 2);
n_frames = size(sim_data_mat(1).traces, 1);
n_reps = size(sim_data_mat(1).traces, 2);
n_odors = size(sim_data_mat(1).traces, 3);
n_durs = size(sim_data_mat(1).traces, 4);
frame_time = sim_data_mat(1).frame_time;
win_end_fr = floor(integration_window./frame_time);

%building odor response matrix as well as mean response matrix.
resp_size_mat = zeros((n_odors.*n_reps), n_cells) + nan;
mean_resp_mat = zeros(n_odors, n_cells) + nan;
for cell_n = 1:n_cells
    for odor_n = 0:(n_odors - 1)
        curr_traces = squeeze(sim_data_mat(cell_n).traces(:, :, (odor_n+1), duration));
        stim_frames = sim_data_mat(cell_n).stim_frs(:, duration);
        resps = nanmean(curr_traces(stim_frames(1):(stim_frames(1) + win_end_fr), :), 1);
        resp_size_mat( ((odor_n.*n_reps + 1):((odor_n + 1).*n_reps) ), cell_n) = resps;  
        mean_resp_mat((odor_n + 1), cell_n) = nanmean(resps);
    end
end


%%
%creating a template-matching decoder. This is to quantify odor separability in the KC representation
test_reps = zeros(n_odors, n_cells) + nan;
template_matrix = zeros(n_odors, n_cells) + nan;
for odor_n = 0:(n_odors - 1)
    curr_resps = resp_size_mat( ((odor_n.*n_reps)+1):((odor_n + 1).*n_reps) , :);
    
    %randomly choosing one repeat to set aside for testing decoder layer performance
    not_nan = 0;
    while not_nan == 0
        picked_rep_n = round(rand(1, 1).*(n_reps - 1) + 1);
        other_reps = curr_resps;
        if nansum(isnan(curr_resps(picked_rep_n, :))) < n_cells./2;
            picked_rep = curr_resps(picked_rep_n, :);
            test_reps( (odor_n + 1), :) = picked_rep;
            other_reps(picked_rep_n, :) = nan;
            not_nan = 1;
        else
            not_nan = 0;
        end
    end
    
    template_matrix( (odor_n + 1), :) = nanmean(other_reps);
    del = isnan(template_matrix);
    del = find(del == 1);
    template_matrix(del) = 0;
    
end

%checking template scores for each of the test trials
for test_trial_n = 1:size(test_reps, 1)
    curr_rep = test_reps(test_trial_n, :);
    del = isnan(curr_rep);
    del = find(del == 1);
    curr_rep(del) = 0;
    corr_mat = corrcoef([curr_rep; template_matrix]');
    corr_vec = corr_mat(1, :);
    keyboard
end
%PICK UP THREAD HERE
%score template match and decoder performance, over-all. 

classifier_output = [];


end