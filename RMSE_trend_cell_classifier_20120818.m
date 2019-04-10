function[test_lists PSTH_mat_full pki] = RMSE_trend_cell_classifier_20120818(dff_data_mat, pk_behav_trial, blink_list, kernel_width, cell_no_fraction, CS_onset_frame, US_onset_frame, frame_time, blink_trials_only)
%RMSE_trend_cell_classifier_20120818.m classifies cells using three
%different criteria - as cells stably encoding time fields. All three
%criteria (or 'tests') use a basic RMSE score calculated by comparing individual trial
%traces with that obtained by averaging trials. 



no_cells = size(dff_data_mat, 2);
no_trials = size(dff_data_mat, 3);

%calculating PSTH mats and pks from pk_behav_trial onwards to end
PSTH_mat_full = nanmean(dff_data_mat(:, :, pk_behav_trial:no_trials), 3);
PSTH_mat = nanmean(dff_data_mat(CS_onset_frame:(US_onset_frame + round(200./frame_time)), :, pk_behav_trial:no_trials ), 3);
[pks pki] = nanmax(PSTH_mat, [], 1);
pki = pki + CS_onset_frame - 1;


%calculating PSTH mat and pks only from blink trials
blink_trials = find(blink_list == 1);
if length(blink_trials)<size(dff_data_mat, 3) || blink_trials_only == 0
    blink_trials = pk_behav_trial:no_trials;
else
end


PSTH_mat_full_blinks = nanmean(dff_data_mat(:, :, blink_trials), 3);
PSTH_mat_blinks = nanmean(dff_data_mat(CS_onset_frame:(US_onset_frame + round(200./frame_time)), :, blink_trials), 3);
[pks_blink pki_blinks] = nanmax(PSTH_mat_blinks, [], 1);
pki_blinks = pki_blinks + CS_onset_frame - 1;


%loops to calculate RMSEs for a small kernel around pk for each cell
RMSE_mat = zeros(no_cells, no_trials);
RMSE_mat_blinks = zeros(no_cells, no_trials);


%kernel_width = 1;
test_lists = zeros(no_cells, 3);
for cell_no = 1:no_cells
    kernel = PSTH_mat_full( (pki(cell_no) - kernel_width):(pki(cell_no) + kernel_width), cell_no); 
    kernel_blinks = PSTH_mat_full_blinks( (pki_blinks(cell_no) - kernel_width):(pki_blinks(cell_no) + kernel_width), cell_no); 
    traces = squeeze(dff_data_mat(:, cell_no, :));
    traces_blinks = squeeze(dff_data_mat(:, cell_no, :));

    RMSE_vec = zeros(1, no_trials);
    RMSE_vec_blinks = zeros(1, no_trials);
    for trial_no = 1:no_trials
        vec = traces((pki(cell_no) - kernel_width):(pki(cell_no) + kernel_width), trial_no);
        RMSE_vec(1, trial_no) = sqrt(mean( (kernel - vec).^2 ) );
        RMSE_vec_blinks(1, trial_no) = sqrt(mean( (kernel_blinks - vec).^2 ) );
    end



    %Tests to look for cells with appropriate RMSEs
    mid = pk_behav_trial - 5;
%     means = [mean(RMSE_vec(1:mid)), mean(RMSE_vec(mid:length(RMSE_vec) ))];
%     stds = [std(RMSE_vec(1:mid)), std(RMSE_vec(mid:length(RMSE_vec) )) ];
%     ses = stds./sqrt(length(RMSE_vec)./2);

    %test1 - Sig Decrease in RMSE after pk behav trial
    [h, p] = ttest2(RMSE_vec(1:mid), RMSE_vec(mid:length(RMSE_vec)) );
    if h == 1
        test_lists(cell_no, 1) = 1;
    else
    end


    RMSE_mat(cell_no, :) = RMSE_vec;
    RMSE_mat_blinks(cell_no, :) = RMSE_vec_blinks;
    clear traces
    clear vec

end

%cell_no_fraction = 0.3;
cell_no_cutoff = round(no_cells.*cell_no_fraction);

%test2 - Taking fraction of cells with lower RMSE means in trials after pk behav trial
means = mean(RMSE_mat(:, pk_behav_trial:no_trials), 2);
means = [means, (1:1:no_cells)'];
means = sortrows(means);
temp = means(1:cell_no_cutoff, 2);
test_lists(temp, 2) = 1;

%test3 - Taking fraction of cells with lower RMSE stds in trials after pk behav trial
stds = std(RMSE_mat(:, pk_behav_trial:no_trials), [], 2);
stds = [stds, (1:1:no_cells)'];
stds = sortrows(stds);
temp = stds(1:cell_no_cutoff, 2);
test_lists(temp, 3) = 1;

clear means
clear stds
        
