function [params_spec1] = setup_params_pairing_expt_elig_trace_Yoshi()
%syntax: function [params_spec1] = setup_params_pairing_expt(paired_odours, led_elec, CS_dur, US_dur, US_dutycyc)
%This function sets up a detailed stimulus specification structure and saves
%it into curr_aq_direc to set up stimulus delivery for a pre, pairing and
%post experiment protocol. Paired odor is the vector of odor numbers to be paired.
%if led_elec is 0, led is paired, if led_elec is 1, elec is paired.

%Reading in last foreground odour to ensure alternate long-pulse odours in each session for reciprocal analysis.
last_paired_od = load('E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\reciprocal_expt_logs\pairing_expt_elig_trace_Yoshi.mat');
last_paired_od = last_paired_od.curr_paired_od;

%step1:Setting up the pre trials
%step1.1 setting up pre-post odours to be delivered with olf1
%creating param specification structure for std KC recording expt. This is for pre-pairing characterisation. 
params_spec1 = set_ADO_params_elig_trace;

params_spec1.reps = 6;
params_spec1.duration = 5;
params_spec1.isi = 52;
params_struc1 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.

%step2: Setting up pairing trial and CS- trial.
%editing param specification structure for only pairing trial. 
params_spec2 = params_spec1;

if last_paired_od == 1
    params_spec2.odours = 4;
    curr_paired_od = 4;
elseif last_paired_od == 4
    params_spec2.odours = 1;
    curr_paired_od = 1;
end
params_spec2.duration = 60;
params_spec2.reps = 1;
params_spec2.randomize = 0;
params_spec2.isi = 250;
params_spec2.led_odours = params_spec2.odours;
params_spec2.stim_dur = 60;
params_spec2.rel_stim_init_delay = 120;
params_spec2.stim_freq = 1;
params_spec2.st_duty_cyc = 5;       %5 percent
params_struc2a = setUpStimuli_modular(params_spec2);


%adding CS- trial
params_spec2.led_odours = [];
params_spec2.odours = last_paired_od;
params_spec2.isi = 120;
params_struc2b = setUpStimuli_modular(params_spec2);


%concatenating pairing trials to params
params = append_params(params_struc1, params_struc2a, 0);  
params = append_params(params, params_struc2b, 0);


%step3: Adding on the post trials, a copy of the pre trials
params = append_params(params, params_struc1, 0);
params_mat = params;

%saving params to current acquisition directory
curr_dir = curr_aq_direc;
if exist([curr_dir, 'params.mat']) == 2
    ovwrite = input('paramfile already exists - overwrite (0-no, 1-yes)?');
    if ovwrite == 1
        save([curr_dir, 'params.mat'], 'params_mat');
        %keeping track of which odor (Hep or IAA) was the foreground odor in this experiment and saving to file. This will be used to alternate the CS+ odor.
        last_long_od = params_spec2.odours;     %This is the long-pulse odour delivered on olf1
        save('E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\reciprocal_expt_logs\pairing_expt_elig_trace_Yoshi.mat', 'curr_paired_od');

    elseif ovwrite == 0
    end
else
    save([curr_dir, 'params.mat'], 'params_mat');
end

disp(['saved detailed stim params structure in ', curr_dir]);



