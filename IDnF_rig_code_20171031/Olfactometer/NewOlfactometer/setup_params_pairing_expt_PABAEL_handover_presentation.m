function [params_spec1] = setup_params_pairing_expt_PABAEL_handover_presentation()
%syntax: function [params_spec1] = setup_params_pairing_expt(paired_odours, led_elec, CS_dur, US_dur, US_dutycyc)
%This function sets up a detailed stimulus specification structure and saves
%it into curr_aq_direc to set up stimulus delivery for a pre, pairing and
%post experiment protocol. Paired odor is the vector of odor numbers to be paired.
%if led_elec is 0, led is paired, if led_elec is 1, elec is paired.

%%Reading in last foreground odour to ensure alternate paired odours in each session for reciprocal analysis.
log_file_path = 'E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\reciprocal_expt_logs\pairing_expt_PABAEL_handover_presentation.mat';
last_paired_od = load(log_file_path);
last_paired_od = last_paired_od.last_paired_od;
odors_to_switch = [3, 10];     %PA and BA on olf1
od_remi = find(odors_to_switch == intersect(last_paired_od, odors_to_switch));
odors_to_switch(od_remi) = [];
od_to_pair = odors_to_switch;
last_paired_od = od_to_pair;
save(log_file_path, 'last_paired_od');

%step1:Setting up the pre trials
%step1.1 setting up PABA together
%creating param specification structure for std KC recording expt. This is for pre-pairing characterisation. 
params_spec1 = set_ADO_params_default;
%since this is a PABAEL expt, odourNames_olf2 is non-default. changing.
del = input('WARNING, PABAEL experiment - swap in BA as odor 3 on olfactometer 2!; Press Enter to conitnue.');
[del, odourNames_olf2]=xlsread('E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2 _PABAEL.xls');
params_spec1.odourNames_olf2 = odourNames_olf2;


%PICK UP THREAD HERE: ensure od_to_pair is the paired odor, fix all the
%other protocol parameter specifications.

keyboard
params_spec1.odours = [3];%PA
params_spec1.odours_olf2 = [2];%BA
params_spec1.reps = 5;
params_spec1.duration = 10;
params_spec1.isi = 60;
params_spec1.duration_olf2 = 10;
params_spec1.rel_stimLatency_olf2 = 10;
params_struc1 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.

%step1.2 setting up PA alone (olf1)
params_spec1 = set_ADO_params_default;

params_spec1.odours = [3];%PA
params_spec1.odours_olf2 = [2];%BA
params_spec1.reps = 5;
params_spec1.duration = 10;
params_spec1.isi = 60;
params_spec1.duration_olf2 = 0;
params_spec1.rel_stimLatency_olf2 = 10;

params_struc2 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.

%step1.3 setting up BA alone (olf2)
params_spec1 = set_ADO_params_default;

params_spec1.odours = [3];%PA
params_spec1.odours_olf2 = [2];%BA
params_spec1.reps = 5;
params_spec1.duration = 0.1000;
params_spec1.isi = 60;
params_spec1.duration_olf2 = 10;
params_spec1.rel_stimLatency_olf2 = 1;

params_struc3 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.

%step1.4 setting up EL alone control (olf1)
params_spec1 = set_ADO_params_default;

params_spec1.odours = [11];%EL
params_spec1.odours_olf2 = [2];%BA
params_spec1.reps = 5;
params_spec1.duration = 10;
params_spec1.isi = 60;
params_spec1.duration_olf2 = 0;
params_spec1.rel_stimLatency_olf2 = 1;

params_struc4 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.


%combining thre results of steps 1.1, 1.2, 1.3 and 1.4
params_struc_pre = append_params(params_struc1, params_struc2,  1);
params_struc_pre = append_params(params_struc_pre,params_struc3, 1);
params_struc_pre = append_params(params_struc_pre, params_struc4, 1);
params = params_struc_pre;           %this is the main parameter structure that will ultimately be saved in curr_aq_direc.


%step2: Setting up the pairing trial
%editing param specification structure for only pairing trial. 
params_spec2 = set_ADO_params_fore_distracter_HepIAA;

 if last_long_od == 5    %case when last expt was with 2-Hep as foreground odour
    params_spec2.odours = 3;
    params_spec2.duration = 60;
    params_spec2.led_odours = 3;
    params_spec2.odours_olf2 = 2;
    params_spec2.duration_olf2 = 60;
    params_spec2.rel_stimLatency = 90;
    params_spec2.stim_dur = 60;
    
else
end

params_struc_pairing = setUpStimuli_modular(params_spec2);         %detailed, trial-by-trial parameter specification structure.

%concatenating pairing trial to params
params = append_params(params, params_struc_pairing, 0);


%step3: Adding on the post trials
params = append_params(params, params_struc_pre, 0);
params_mat = params;

%saving params to current acquisition directory
curr_dir = curr_aq_direc;
if exist([curr_dir, 'params.mat']) == 2
    ovwrite = input('paramfile already exists - overwrite (0-no, 1-yes)?');
    if ovwrite == 1
        save([curr_dir, 'params.mat'], 'params_mat');
        %keeping track of which odor (Hep or IAA) was the foreground odor in this experiment and saving to file. This will be used to alternate the CS+ odor.
%         last_long_od = params_spec2.odours;     %This is the long-pulse odour delivered on olf1
%         save('E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\reciprocal_expt_logs\pairing_expt_integrator_HepIAA_lastlongod.mat', 'last_long_od');

    elseif ovwrite == 0
    end
else
    save([curr_dir, 'params.mat'], 'params_mat');
end

disp(['saved detailed stim params structure in ', curr_dir]);



