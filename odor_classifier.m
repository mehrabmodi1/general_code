function [classifier_output] = odor_classifier(sim_data_mat, duration, integration_window)

n_cells = size(sim_data_mat, 2);
n_frames = size(sim_data_mat(1).traces, 1);
n_reps = size(sim_data_mat(1).traces, 2);
n_odors = size(sim_data_mat(1).traces, 3);
n_durs = size(sim_data_mat(1).traces, 4);
frame_time = sim_data_mat(1).frame_time;
win_end_fr = floor(integration_window./frame_time);

%building odor response matrix as well as mean respnse matrix to run a PCA on.
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

%Running PCA
%[pc_weights, scores, latent, t_squared, explained] = pca(resp_size_mat, 'algorithm','als');
[pc_weights, scores1, latent, t_squared, explained] = pca(mean_resp_mat);

% figure(1)
% subplot(2, 1, 1)
% imagesc(pc_weights)
% title('pc weight matrix')
% ylabel('cell number')
% xlabel('PC number')
% subplot(2, 1, 2)
% imagesc(resp_size_mat)
% title('dF/F response size matrix')
% ylabel('sorted trial number')
% xlabel('cell number')

%re-calculating scores ignoring nans in resp_size_mat. doing this because
%if matlab comes across a nan for a single cell, it removes an entire row
%from the scores matrix
n_pcs = size(pc_weights, 2);
scores = zeros(size(scores1, 1), size(scores1, 2));
for trial_n = 1:size(resp_size_mat, 1)
    for pc_n = 1:n_pcs
        curr_weights = pc_weights(:, pc_n);
        curr_resps = resp_size_mat(trial_n, :)';
        scores(trial_n, pc_n) = nanmean(curr_weights.*curr_resps);
    end
end
scores = scores.*100;     %rescaling to match scores matrix calculated by matlab's pca function

%%
%Finding centroid positions for each odor cloud in PC space, also measuring
%cloud size as mean distance from centroid for each cloud. Then measuring
%distances between pairs of centroids normalised to the sums of their cloud
%widths.
centroids = zeros(n_odors, n_pcs) + nan;
for odor_n = 0:(n_odors - 1)
    curr_trs = (odor_n.*n_reps + 1):( (odor_n+ 1).*n_reps);      %list of trials for current odor
    curr_scores = scores(curr_trs, :);                           %PC scores for current odor trials
    centroids((odor_n + 1), :) = nanmean(curr_scores, 1);        %centroids for repeats in current PC space. should be similar to PC scores for mean dF/F trace.
    
    %calculating distance of each trial point from the centroid
    curr_scores = [centroids( (odor_n + 1), :); curr_scores];
    curr_dists = pdist(curr_scores);
    curr_dists = curr_dists(1, 1:n_reps);
    cloud_sds( (odor_n + 1), 1 ) = nanmean(curr_dists);                     %mean distance of points from centroid is SD of cloud around it's mean (centroid).
end

%calculating pairwise distances between centroids
cent_dists = pdist(centroids);

%creating pairwise sum(cloud SDs) vector to normalise centroid distances
pairwise_sd_sums = [];
for odor_n = 1:(n_odors-1)
    for odor_ni = (odor_n + 1):n_odors
        pairwise_sd_sums = [pairwise_sd_sums; (cloud_sds(odor_n) + cloud_sds(odor_ni))];
    end
end


cent_dists = cent_dists./pairwise_sd_sums';      %computing z-scored centroid distances
cent_dists = squareform(cent_dists);

%calculating a measure of cloud separability (mean z-score)
del = eye(size(cent_dists, 1));
del = find(del == 1);
cent_dists(del) = nan;
mean_dist = nanmean(nanmean(cent_dists));
fail_frac = length(find(cent_dists < 1))./(size(cent_dists, 1).^2);        

% figure(1)
% sub_scores(:, :) = centroids(:, 1:3)'; 
% cloud_plotter(1, sub_scores, cloud_sds, .3);
% title(['n PCs = ' int2str(n_pcs)])
%PICK UP THREAD HERE
%plot PCA plots and do a reality check. 

classifier_output = [mean_dist, fail_frac];

end