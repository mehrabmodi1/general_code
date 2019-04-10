function [classifier_output] = odor_classifier(sim_data_mat, duration, integration_window, sp_vec)

n_cells = size(sim_data_mat, 2);
n_frames = size(sim_data_mat(1).traces, 1);
n_reps = size(sim_data_mat(1).traces, 2);
n_odors = size(sim_data_mat(1).traces, 3);
n_durs = size(sim_data_mat(1).traces, 4);
frame_time = sim_data_mat(1).frame_time;
win_end_fr = floor(integration_window./frame_time);

%building odor response matrix as well as mean response matrix.
resp_size_mat = zeros((n_odors.*n_reps), n_cells) + nan;
resp_size_mat_thresh = resp_size_mat;
for cell_n = 1:n_cells
    for odor_n = 0:(n_odors - 1)
        curr_traces = squeeze(sim_data_mat(cell_n).traces(:, :, (odor_n+1), duration));
        stim_frames = sim_data_mat(cell_n).stim_frs(:, duration);
        resps = nanmean(curr_traces(stim_frames(1):(stim_frames(1) + win_end_fr), :), 1);
        base_sd = nanstd(curr_traces(1:stim_frames(1), :), [], 1);
        resp_size_mat( ((odor_n.*n_reps + 1):((odor_n + 1).*n_reps) ), cell_n) = resps;  
        
        %creating a thresholded mean response matrix
        resps_sub = resps - (base_sd);
        resps_sub = sign(resps_sub);
        resps_sub(resps_sub < 1) = 0;
        resp_size_mat_thresh(((odor_n.*n_reps + 1):((odor_n + 1).*n_reps) ), cell_n) = resps.*resps_sub;
        %NOTE: tried using the thresholded response matrix as input to the
        %template matching decoder and performance was much worse. Unclear
        %why this was the case.
    
    end
end


%%
%creating a template-matching decoder. This is to quantify odor separability in the KC representation
test_reps = zeros(n_odors, n_cells) + nan;
template_matrix = zeros(n_odors, n_cells) + nan;

for rep_n = 1:n_reps
    for odor_n = 0:(n_odors - 1)
        curr_resps = resp_size_mat( ((odor_n.*n_reps)+1):((odor_n + 1).*n_reps) , :);
        test_reps((odor_n + 1), :) = curr_resps(rep_n, :);
        
        curr_resps(rep_n, :) = [];
        template_matrix( (odor_n + 1), :) = nanmean(curr_resps);
        del = isnan(template_matrix);
        del = find(del == 1);
        template_matrix(del) = 0;

    end

    %checking template scores for each of the test trials
    preds_mat = zeros(n_odors, n_reps);
    for test_trial_n = 1:size(test_reps, 1)
        curr_rep = test_reps(test_trial_n, :);
        del = isnan(curr_rep);
        del = find(del == 1);
        curr_rep(del) = 0;
        corr_mat = corrcoef([curr_rep; template_matrix]');
        corr_vec = corr_mat(1, :);          %correlations with test trial are only the first row or col of corr mat
        corr_vec(1) = [];                   %getting rid of correlation with self
        [del, pred_odor] = nanmax(corr_vec);
        pred_mat(test_trial_n, rep_n) = pred_odor;

        
    end

end

%scoring performance
%creating confusion matrix
confusion_mat = zeros(n_odors, n_odors);
for odor_n = 1:n_odors
    for rep_n = 1:n_reps
        curr_guess = pred_mat(odor_n, rep_n);
        confusion_mat(odor_n, curr_guess) = confusion_mat(odor_n, curr_guess) + 1; 
    end
end
confusion_mat = confusion_mat./n_reps;
correct_mat = confusion_mat.*eye(n_odors, n_odors);
frac_correct = sum(sum(correct_mat))./n_odors;

% figure(1)
% imagesc(confusion_mat)
% xlabel('odor number presented')
% ylabel('odor number classified')

odor_fail_fracs = sum(correct_mat, 1);                  %fraction of correct guesses for each odor
c = corrcoef(odor_fail_fracs, sp_vec);
c = c(1, 2);
% figure(2)
% plot(odor_fail_fracs, sp_vec, '.')
% title(['corrcoef = ' num2str(c, 3)])
% xlabel('fraction of correct classifications')
% ylabel('odor response sparseness')

classifier_output = [{confusion_mat}, {frac_correct}];
end