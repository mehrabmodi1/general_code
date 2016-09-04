function [dff_all_PCs, areas_all_PCs] = PCA_calculator(resp_areas, nPCs, dff_data_mat, odor_list);
%This function calculates principal components using the response areas on
%all trials, it then re-calculates dff data traces for all trials by
%re-weighting cells by PC coefficients and averaging across these to get
%each PC response time-series. Its output dff_all_PCs is a n_frames by n_PCs by
%n_trials by n_odors matrix.
%syntax: [dff_all_PCs] = PCA_calculator(resp_areas, nPCs, dff_data_mat, odor_list);
%Mehrab Modi, 20150803

n_frames = size(dff_data_mat, 1);
n_cells = size(dff_data_mat, 2);
n_trials = size(dff_data_mat, 3);
n_odors = length(odor_list);

dff_reweighted = zeros(n_frames, n_cells, n_trials, n_odors) + nan;
dff_all_PCs = zeros(n_frames, nPCs, n_trials, n_odors) + nan;

X = resp_areas';        %calculating PCs on only response areas, not full time series
[coeff,score,latent,tsquared,explained] = pca(X);
areas_allPCs = zeros(nPCs, n_trials);

for PC_n = 1:nPCs
    
    coeffs = coeff(:, PC_n);
    
    
    coeff_mat = repmat(coeffs', [n_frames, 1, n_trials, 8]);

    %recalculating dff_data_mat for the current PC
    dff_reweighted = dff_data_mat.*coeff_mat;
    dff_PCn = nanmean(dff_reweighted, 2);                     %averaging reweighted dff traces to get PC projection for current PC
    dff_reweighted = zeros(n_frames, n_cells, n_trials, n_odors) + nan;
    dff_all_PCs(:, PC_n, :, :) = dff_PCn(:, :, :, odor_list);                  %storing projection for current PC
    
    %recalculating resp_areas for the current PC
    coeff_mat = repmat(coeffs, [1, n_trials]);
    resp_areas_reweighted = resp_areas.*coeff_mat;
    areas_PCn = nanmean(resp_areas_reweighted, 1);
    areas_all_PCs(PC_n, :) = areas_PCn;
    
end
