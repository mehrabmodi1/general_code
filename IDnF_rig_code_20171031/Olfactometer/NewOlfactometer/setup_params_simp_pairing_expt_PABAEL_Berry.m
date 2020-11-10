function [params_spec1] = setup_params_simp_pairing_expt_PABAEL_Berry(LED_on)
%syntax: [params_spec1] = setup_params_simp_pairing_expt_PABAEL(LED_on)
%This function sets up a detailed stimulus specification structure and saves
%it into curr_aq_direc to set up stimulus delivery for a protocol matching
%the Berry et al. 2018 paper.

%%Reading in last paired odour to ensure alternate paired odours in each session for reciprocal analysis.
log_file_path = 'E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\reciprocal_expt_logs\pairing_expt_PABAEL_simple_presentation.mat';
last_paired_od = load(log_file_path);
last_paired_od = last_paired_od.last_paired_od;
odors_to_switch = [10, 3];     %PA and BA on olf2
od_remi = find(odors_to_switch == intersect(last_paired_od, odors_to_switch));
odors_to_switch(od_remi) = [];
CSplus_od = odors_to_switch;
CSminus_od = last_paired_od;
last_paired_od = CSplus_od;     %swapping to log for next expt.
save(log_file_path, 'last_paired_od');

%step1:Setting up the pre trials
%step1.1 setting up pre trials - no imaging
%creating param specification structure for std KC recording expt. This is for pre-pairing characterisation. 
params_spec1 = set_ADO_params_default;
%since this is a PABAEL expt, odourNames_olf2 is non-default. changing.
del = input('WARNING, PABAEL experiment - swap in BA as odor 3 on olfactometer 2!; Press Enter to conitnue.');
[del, odourNames_olf2]=xlsread('E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2 _PABAEL.xls');
params_spec1.odourNames_olf2 = odourNames_olf2;
params_spec1.odours = [3, 10, 11];
params_spec1.duration = 5;
params_spec1.reps = 1;
params_spec1.isi = 45;
params_spec1.odours_olf2 = [];
params_spec1.trigger_scan = 0;

params_struc1 = setUpStimuli_modular(params_spec1);

%step1.2 setting up pre trials - with imaging
params_spec1.trigger_scan = 1;
params_struc1_1 = setUpStimuli_modular(params_spec1);
params_struc_pre = append_params(params_struc1, params_struc1_1, 0);

%step2.1: Setting up the pairing trial
%editing param specification structure for only pairing trial. 
params_spec2 = set_ADO_params_fore_distracter_HepIAA;

params_spec2.odours = CSplus_od;
params_spec2.duration = 5;
if LED_on == 0
    params_spec2.led_odours = [];
else
    params_spec2.led_odours = CSplus_od;
end
params_spec2.n_od_pulses = 1;
params_spec2.stim_dur = 3;
params_spec2.stim_freq = 1;
params_spec2.st_duty_cyc = 100;
params_spec2.rel_stim_init_delay = 2;
params_spec2.isi = 45;     
params_spec2.trigger_scan = 0;
params_spec2.odours_olf2 = [];

params_struc_pairing = setUpStimuli_modular(params_spec2);         %detailed, trial-by-trial parameter specification structure.

%Step2.2 Setting up the CS- trial
params_spec3 = params_spec1;
params_spec3.odours = [CSminus_od, 11];
params_spec3.reps = 1;
params_spec3.trigger_scan = 0;

params_struc_CSminus = setUpStimuli_modular(params_spec3);         %detailed, trial-by-trial parameter specification structure.

%looping the pairing and CS- trials thrice
params_struc = params_struc_pre;
for rep_n = 1:3
    params_struc = append_params(params_struc, params_struc_pairing, 0);
    params_struc = append_params(params_struc, params_struc_CSminus, 0);
end

%step3: Adding on the post1 trials, first the imaged trials and then
%un-imaged trials
params_struc = append_params(params_struc, params_struc1_1, 0);
params_struc = append_params(params_struc, params_struc1, 0);

%step4: Adding on un-paired LED stim trial for CS+ odor
params_spec4 = params_spec2;
params_spec4.rel_stim_init_delay = 22.5;
params_struc_unpairing = setUpStimuli_modular(params_spec4);

%looping the un-pairing and CS- trials thrice
for rep_n = 1:3
    params_struc = append_params(params_struc, params_struc_unpairing, 0);
    params_struc = append_params(params_struc, params_struc_CSminus, 0);
end

%Step5: Adding on the post2 trials, first un-imaged trials and then imaged
%trials
params_struc = append_params(params_struc, params_struc1, 0);
params_struc = append_params(params_struc, params_struc1_1, 0);

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
