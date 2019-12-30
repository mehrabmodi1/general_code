function [params_spec1] = setup_params_KC_PABAEL_handover_presentation_orig()
%syntax: [params_spec1] = setup_params_pairing_expt_PABAEL_handover_presentation()
%This function sets up a detailed stimulus specification structure and saves
%it into curr_aq_direc to set up stimulus delivery for a pre, pairing and
%post experiment protocol, with extra trials to hand off a CS+ to CS- odor pulses and vice
%versa. 

%%Reading in last olf2 odour to ensure alternate paired odours in each session for reciprocal analysis.
log_file_path = 'E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\reciprocal_expt_logs\KC_PABAEL_handover_presentation.mat';
% last_paired_od = load(log_file_path);
% last_paired_od = last_paired_od.last_paired_od;
last_paired_od = 1;
odors_to_switch = [1, 3];     %PA and BA on olf2
od_remi = find(odors_to_switch == intersect(last_paired_od, odors_to_switch));
odors_to_switch(od_remi) = [];
od_to_pair = odors_to_switch;
last_paired_od = od_to_pair;
save(log_file_path, 'last_paired_od');

%step1:Setting up the odor trials
%step1.1 setting up PABA together
params_spec1 = set_ADO_params_default;
%since this is a PABAEL expt, odourNames_olf2 is non-default. changing.
del = input('WARNING, PABAEL experiment - swap in BA as odor 3 on olfactometer 2!; Press Enter to conitnue.');
[del, odourNames_olf2]=xlsread('E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2 _PABAEL.xls');
params_spec1.odourNames_olf2 = odourNames_olf2;

%the paired odor is presented on olf2.
if last_paired_od == 1 %PA on olf2
    params_spec1.odours = 3;        %PA on olf1
    params_spec1.odours_olf2 = 3;   %BA on olf2
elseif last_paired_od == 3 %BA on olf2
    params_spec1.odours = 10;       %BA on olf1
    params_spec1.odours_olf2 = 1;   %PA on olf2
else
end    

%step1.1 setting up olf1 alone trials
od_olf2 = params_spec1.odours_olf2;
params_spec1.odours_olf2 = [];%BA
params_spec1.reps = 5;
params_spec1.duration = 10;
params_spec1.isi = 110;
params_spec1.post_od_scan_dur = 70;
params_struc0 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.

%step 1.2 setting up olf2 alone trials
params_spec1.odours_olf2 = od_olf2;%BA
params_spec1.duration = 0.1;
params_spec1.duration_olf2 = 10;
params_spec1.rel_stimLatency_olf2 = 0;
params_struc1 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.

%This contains olf1 alone and olf2 alone trials
params_struc = append_params(params_struc0, params_struc1, 1);  %combining and randomising explicit param spec structures


%step 1.3 adding handover olf1-olf2 trials (CS- to CS+ odor trials)   
params_spec1.reps = 5;
params_spec1.duration = 10;
params_spec1.isi = 60;
params_spec1.duration_olf2 = 10;
params_spec1.rel_stimLatency_olf2 = 10;
params_spec1.post_od_scan_dur = 10;
params_struc2 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.

params_struc = append_params(params_struc, params_struc2, 1);  %combining and randomising explicit param spec structures

%step 1.4 adding handover olf1-olf2 trials (CS+ to CS- odor trials)   
params_spec2 = params_spec1;
%the paired odor is presented on olf2.
if last_paired_od == 1 %PA on olf2
    params_spec2.odours = 10;%BA
    params_spec2.odours_olf2 = 1;%PA
elseif last_paired_od == 3 %BA on olf2
    params_spec2.odours = 3;%PA
    params_spec2.odours_olf2 = 3;%BA
else
end    
params_struc3= setUpStimuli_modular(params_spec2);         %detailed, trial-by-trial parameter specification structure for olf1 odour.
params_struc = append_params(params_struc, params_struc3, 1);  %combining and randomising explicit param spec structures


%step1.5 setting up EL alone control (olf1)
params_spec1 = set_ADO_params_default;

params_spec1.odours = [11];%EL
params_spec1.odours_olf2 = [];
params_spec1.reps = 5;
params_spec1.duration = 10;
params_spec1.isi = 110;
params_spec1.post_od_scan_dur = 70;
params_struc3 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.

params_struc = append_params(params_struc, params_struc3, 1);  %combining and randomising explicit param spec structures
params_mat = params_struc;

%saving params to current acquisition directory
curr_dir = curr_aq_direc;
if exist([curr_dir, 'params.mat']) == 2
    ovwrite = input('paramfile already exists - overwrite (0-no, 1-yes)?');
    if ovwrite == 1
        save([curr_dir, 'params.mat'], 'params_mat');
        
    elseif ovwrite == 0
    end
else
    save([curr_dir, 'params.mat'], 'params_mat');
end

disp(['saved detailed stim params structure in ', curr_dir]);


