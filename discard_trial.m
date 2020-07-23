function [] = discard_trial(direc, trial_n)
%syntax: function [] = discard_trial(direc, trial_n)
%This function allows the user to manually identify trials with issues
%(movement etc)and discard them in the saved_traces .mat file

extracted_raw_data_mat = load([direc, 'extracted_raw_data_mat.mat']);
raw_data_mat = extracted_raw_data_mat.raw_data_mat;
raw_data_mat_orig = raw_data_mat;

if isempty(trial_n) == 0
    raw_data_mat(:, :, trial_n) = nan;
    
elseif isempty(trial_n) == 1    %case where user wants to revert to original, full raw_data_mat
    extracted_raw_data_mat = load([direc, 'extracted_raw_data_mat_orig.mat']);
    raw_data_mat = extracted_raw_data_mat.raw_data_mat_orig;
    raw_data_mat_orig = raw_data_mat;
else
end
    
%saving
save([direc, 'extracted_raw_data_mat.mat'], 'raw_data_mat');
save([direc, 'extracted_raw_data_mat_orig.mat'], 'raw_data_mat_orig');
