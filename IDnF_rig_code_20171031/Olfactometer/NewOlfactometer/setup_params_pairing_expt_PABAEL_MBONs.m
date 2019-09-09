function [params_spec1] = setup_params_pairing_expt_PABAEL_MBONs()
%syntax: function [params_spec1] = setup_params_pairing_expt(paired_odours, led_elec, CS_dur, US_dur, US_dutycyc)
%This function sets up a detailed stimulus specification structure and saves
%it into curr_aq_direc to set up stimulus delivery for a pre, pairing and
%post experiment protocol. Paired odor is the vector of odor numbers to be paired.
%if led_elec is 0, led is paired, if led_elec is 1, elec is paired.

% %Reading in last foreground odour to ensure alternate long-pulse odours in each session for reciprocal analysis.
% last_long_od = load('E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\reciprocal_expt_logs\pairing_expt_PABA_EL_lastlongod.mat');
% last_long_od = last_long_od.last_long_od;
last_long_od = 3;

%step1:Setting up the pre trials
%step1.1 setting up PABA together
%creating param specification structure for std KC recording expt. This is for pre-pairing characterisation. 
params_spec1 = set_ADO_params_default;

params_spec1.odours = [3, 10, 11];%PA
params_spec1.odours_olf2 = [];%BA
params_spec1.reps = 7;
params_spec1.duration = 10;
params_spec1.isi = 60;
params_spec1.duration_olf2 = 10;
params_spec1.rel_stimLatency_olf2 = 10;
params_spec1.randomize = 1;
params_struc1 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.
params = params_struc1;

%step2: Setting up the pairing trial
%editing param specification structure for only pairing trial. 
params_spec2 = params_spec1;

if last_long_od == 3    %case when last expt was with 2-Hep as foreground odour
    curr_long_od = 10;
elseif last_long_od == 10
    curr_long_od = 3;
end

params_spec2.reps = 1;
params_spec2.isi = 130;  %10 + 60 + 10 + 50
params_spec2.odours = curr_long_od;
params_spec2.duration = 60;
params_spec2.led_odours = curr_long_od;
params_spec2.odours_olf2 = [];
params_spec2.led_odours = curr_long_od;
params_spec2.stim_dur = 60;
params_spec2.st_duty_cyc = 0.5;
params_spec2.rel_stim_init_delay = 5;

params_struc_pairing = setUpStimuli_modular(params_spec2);         %detailed, trial-by-trial parameter specification structure.

%concatenating pairing trial to params
params = append_params(params, params_struc_pairing, 0);

%adding on CS- presentation
params_spec2.odours = last_long_od;
params_spec2.led_odours = [];
params_struc_CSm = setUpStimuli_modular(params_spec2);         %detailed, trial-by-trial parameter specification structure.
params = append_params(params, params_struc_CSm, 0);

%step3: Adding on the post trials

params = append_params(params, params_struc1, 0);
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
    last_long_od = curr_long_od;
    save('E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\reciprocal_expt_logs\pairing_expt_PABAEL_lastlongod.mat', 'last_long_od');
end

disp(['saved detailed stim params structure in ', curr_dir]);



