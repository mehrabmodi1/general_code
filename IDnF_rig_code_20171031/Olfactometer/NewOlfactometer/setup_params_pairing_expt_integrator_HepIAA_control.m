function [params_spec1] = setup_params_pairing_expt_integrator_HepIAA_control(US_dutycyc, LED_power)
%syntax: function [params_spec1] = setup_params_pairing_expt(paired_odours, led_elec, CS_dur, US_dur, US_dutycyc)
%This function sets up a detailed stimulus specification structure and saves
%it into curr_aq_direc to set up stimulus delivery for a pre, pairing and
%post experiment protocol. Paired odor is the vector of odor numbers to be paired.
%if led_elec is 0, led is paired, if led_elec is 1, elec is paired.

%Reading in last foreground odour to ensure alternate long-pulse odours in each session for reciprocal analysis.
last_long_od = load('E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\reciprocal_expt_logs\pairing_expt_integrator_HepIAA_lastlongod.mat');
last_long_od = last_long_od.last_long_od;

%step1:Setting up the pre trials
%step1.1 setting up pre-post odours to be delivered with olf1
%creating param specification structure for std KC recording expt. This is for pre-pairing characterisation. 
params_spec1 = set_ADO_params_HepIAA;
if last_long_od == 9
    params_spec1.odours = [5, 11];
elseif last_long_od == 5
    params_spec1.odours = [9, 11];
else
end

params_spec1.reps = 5;
params_spec1.duration = 10;
params_spec1.isi = 52;
params_spec1.duration_olf2 = 0;
params_struc1 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.

%step1.2 setting up pre-post odor to be delivered with olf2 (the distr odor)
if last_long_od == 5
    params_spec1.odours_olf2 = 3;
elseif last_long_od == 9
    params_spec1.odours_olf2 = 2;
else
end
params_spec1.odours = last_long_od;
params_spec1.duration_olf2 = 10;
params_spec1.duration = 0.05;
params_struc2 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.

%combining thre results of steps 1.1 and 1.2
params_struc_pre = append_params(params_struc1, params_struc2, 1);
params = params_struc_pre;           %this is the main parameter structure that will ultimately be saved in curr_aq_direc.


%step2: Setting up the pairing trial
%editing param specification structure for only pairing trial. 
params_spec2 = set_ADO_params_fore_distracter_HepIAA_control;

if last_long_od == 5    %case when last expt was with 2-Hep as foreground odour
    params_spec2.odours = 9;
    params_spec2.led_odours = 9;
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
        last_long_od = params_spec2.odours;     %This is the long-pulse odour delivered on olf1
        save('E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\reciprocal_expt_logs\pairing_expt_integrator_HepIAA_lastlongod.mat', 'last_long_od');

    elseif ovwrite == 0
    end
else
    save([curr_dir, 'params.mat'], 'params_mat');
end

disp(['saved detailed stim params structure in ', curr_dir]);



