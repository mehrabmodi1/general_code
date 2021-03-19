function [params_spec1] = setup_params_handover_KC_expt_PABAEL_Berry(n_reps, EL_type, int_pulse_int)
%syntax: [params_spec1] = setup_params_simp_pairing_expt_PABAEL(n_reps, EL_type, int_pulse_int)
%This function sets up a detailed stimulus specification structure and saves
%it into curr_aq_direc to set up stimulus delivery for PA-BA-EL sinple and handover trials
%EL_type controls EL in handover trials. 0 - no EL in handover trials; 1 - EL 
%is the first odor in handover trials; 2 - EL is the second odor in handover trials.

%%Reading in last paired odour to ensure alternate paired odours in each session for reciprocal analysis.
log_file_path = 'E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\reciprocal_expt_logs\just_pres_expt_PABAEL_handover_Berry.mat';
last_paired_od = load(log_file_path);
last_paired_od = last_paired_od.last_paired_od;
%last_paired_od = 4;
odors_to_switch = [9, 4];     %PA and BA on olf1 NOTE:valve3 onlof1 stopped working, so nout mounting PA in valve4.
od_remi = find(odors_to_switch == intersect(last_paired_od, odors_to_switch));
odors_to_switch(od_remi) = [];
CSplus_od = odors_to_switch;    %on olf1
CSminus_od = last_paired_od;    %on olf1
last_paired_od = CSplus_od;     %swapping to log for next expt.
save(log_file_path, 'last_paired_od');

od_lookup_table = [4, 1; 9, 3; 11, 4];     %table of corress od numbers on olf1 and olf2.
plusi = find(od_lookup_table(:, 1) == CSplus_od);
CS_plus_od_olf2 = od_lookup_table(plusi, 2);
minusi = find(od_lookup_table(:, 1) == CSminus_od);
CS_minus_od_olf2 = od_lookup_table(minusi, 2);

%step1:Setting up the pre trials
%step1.1 setting up pre trials - no imaging

%step1.1.1 setting up plus-minus handover trial 
params_spec1 = set_ADO_params_default;
%since this is a PABAEL expt, odourNames_olf2 is non-default. changing.
del = input('WARNING, PABAEL experiment - swap in BA as odor 3 on olfactometer 2!; Press Enter to conitnue.');
[del, odourNames_olf2]=xlsread('E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2 _PABAEL.xls');
params_spec1.odourNames_olf2 = odourNames_olf2;
if EL_type == 0
    params_spec1.odours = CSplus_od;
    params_spec1.odours_olf2 = CS_minus_od_olf2;
elseif EL_type == 1
    params_spec1.odours = 11;
    params_spec1.odours_olf2 = CS_minus_od_olf2;
elseif EL_type == 2
    params_spec1.odours = CSplus_od;
    params_spec1.odours_olf2 = 4;
else
end 
params_spec1.duration = 5;
params_spec1.reps = n_reps;
params_spec1.isi = 45 + int_pulse_int;
params_spec1.rel_stimLatency_olf2 = params_spec1.duration + int_pulse_int;
params_spec1.duration_olf2 = 5;
params_spec1.trigger_scan = 1;

params_struc1 = setUpStimuli_modular(params_spec1);

%step1.1.2 setting up minus-plus handover trial
if EL_type == 0
    params_spec1.odours = CSminus_od;
    params_spec1.odours_olf2 = CS_plus_od_olf2;
elseif EL_type == 1
    params_spec1.odours = 11;
    params_spec1.odours_olf2 = CS_plus_od_olf2;
elseif EL_type == 2
    params_spec1.odours = CSminus_od;
    params_spec1.odours_olf2 = 4;
else
end
params_struc1_1 = setUpStimuli_modular(params_spec1);

params_struc = append_params(params_struc1, params_struc1_1, 1);

%step1.1.3 setting up simple pulse trials
olf2_od_vec = [1, 3, 4];
for olf2_od_n = 1:3
    params_spec1.odours = 11;
    params_spec1.duration = 0.05;
    params_spec1.odours_olf2 = olf2_od_vec(olf2_od_n);
    params_spec1.rel_stimLatency_olf2 = 0;
    params_struc1_2 = setUpStimuli_modular(params_spec1);

    params_struc = append_params(params_struc, params_struc1_2, 1);
end
params = params_spec1;
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

