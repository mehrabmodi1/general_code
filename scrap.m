clear all
close all

direc = 'C:\Data\Data\Analysed_data\Manual_ROI_results\20200716\fly1_G2_PABAEL_handover_noLED_yesdendrites';
stim_mat = load([direc, '\stim_mat.mat']);
stim_mat = stim_mat.stim_mat;

flip1_count = 0;
flip2_count = 0;
for trial_n = 1:size(stim_mat, 2)
    curr_od2 = stim_mat(trial_n).odours_olf2;
    
    if curr_od2 == 1
        stim_mat(trial_n).odours_olf2 = 3;
        flip1_count = flip1_count + 1
    elseif curr_od2 == 3
        stim_mat(trial_n).odours_olf2 = 1;
        flip2_count = flip2_count + 1
    else
    end
    
    
end

save([direc, '\stim_mat.mat'], 'stim_mat');