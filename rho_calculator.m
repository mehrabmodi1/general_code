function [rho pval] = rho_calculator(dff_data_mat, cell_list, t_list1, t_list2, pre_clip, post_clip);
% rho_calculator syntax: [rho] = rho_calculator(dff_data_mat, t_list1, t_list2, cell_list)
%rho_calculator calculates Spearman's Rho Rank correlations for the ranks
%of cells from cell_list, segragating trials into two groups - t_list1 and
%t_list2. kernel_mat is calculated by averaging all the trials in dff_data_mat.
%The rho value compares the rankings for t_list1 and t_list2 with
%each other. If un-specified, t_list1 and t_list2 default to lists of alternate
%trials in the dataset. pre_clip and post_clip are set by default to
%stimulus and trace times, but can be set to any set of frames - if need
%be. they define the set of frames to be considered for identifying time of
%response.


no_frames = size(dff_data_mat, 1);
no_cells = size(dff_data_mat, 2);
no_trials = size(dff_data_mat, 3);
mid_trial = floor(no_trials./2);


if nargin == 1
    cell_list = 1:size(dff_data_mat, 2);
    %generating two lists of alternate trials
    if rem((no_trials - mid_trial), 2) == 0
        tlist1 = mid_trial:2:no_trials;
        tlist2 = (mid_trial+1):2:(no_trials-1);
    elseif rem((no_trials - mid_trial), 2) == 1
        tlist1 = mid_trial:2:(no_trials-1);
        tlist2 = (mid_trial + 1):2:no_trials;
    end
elseif nargin < 3 | isempty(t_list1) == 1 | isempty(t_list2) == 1
    %generating two lists of alternate trials
    if rem((no_trials - mid_trial), 2) == 0
        tlist1 = mid_trial:2:no_trials;
        tlist2 = (mid_trial+1):2:(no_trials-1);
    elseif rem((no_trials - mid_trial), 2) == 1
        tlist1 = mid_trial:2:(no_trials-1);
        tlist2 = (mid_trial + 1):2:no_trials;
    end
  
end

%calculating two sets of averaged traces for all cells - using two sets of trials
PSTH_mat1 = mean(dff_data_mat(pre_clip:post_clip, cell_list, tlist1), 3);
PSTH_mat2 = mean(dff_data_mat(pre_clip:post_clip, cell_list, tlist2), 3);

[C1 pks1] = max(PSTH_mat1, [], 1);
[C2 pks2] = max(PSTH_mat2, [], 1);

[rho, pval] = corr(pks1', pks2', 'type', 'Spearman', 'tail', 'right' );


