function [params_struc1] = setup_params_trace_pairing_expt(trace_interval, n_pairings, no_LED_control)
%syntax: [params_spec1] = setup_params_pairing_expt_PABAEL_handover_presentation()
%This function sets up a detailed stimulus specification structure and saves
%it into curr_aq_direc to set up stimulus delivery for a pre, pairing and
%post experiment protocol, with extra trials to hand off a CS+ to CS- odor pulses and vice
%versa. 

%%Reading in last paired odour to ensure alternate paired odours in each session for reciprocal analysis.
log_file_path = 'E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\reciprocal_expt_logs\trace_pairing_expt.mat';
last_paired_od = load(log_file_path);
last_paired_od = last_paired_od.last_paired_od;
odors_to_switch = [3, 11];     %PA and EL on olf1
od_remi = find(odors_to_switch == intersect(last_paired_od, odors_to_switch));
odors_to_switch(od_remi) = [];
od_to_pair = odors_to_switch;
last_paired_od = od_to_pair;
save(log_file_path, 'last_paired_od');

%step1:Setting up the pre PA, EL trials
%creating param specification structure for std KC recording expt. This is for pre-pairing characterisation. 
params_spec1 = set_ADO_params_default;
params_spec1.odours_olf2 = [];
params_spec1.duration = 10;
params_spec1.isi = 60;
params_spec1.reps = 5;
params_spec1.odours = [3, 11, 9];    %PA, EL and IAA on olf1
params_spec1.n_od_pulses = 1;

params_struc1 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.
params_struc_pre = params_struc1;

%step2.1: Setting up the pairing trial
%editing param specification structure for only pairing trial. 
params_spec2 = set_ADO_params_fore_distracter_HepIAA;
params_spec2.odours = od_to_pair;
params_spec2.n_od_pulses = 1;
params_spec2.duration = 10;
params_spec2.isi = (params_spec2.stim_dur + trace_interval) + 120;
params_spec2.stimLatency = 10;
params_spec2.post_od_scan_dur = 10;
params_spec2.odours_olf2 =[];
if no_LED_control == 1
    params_spec2.stim_dur = [];
    params_spec2.rel_stim_init_delay = [];      %Yoshi defines ISI as difference between odor and LED onsets
    params_spec2.led_odours = [];
else
    params_spec2.stim_dur = 10;
    params_spec2.rel_stim_init_delay = trace_interval;      %Yoshi defines ISI as difference between odor and LED onsets
    params_spec2.led_odours = od_to_pair;
end
params_spec2.stim_freq = 1;
params_spec2.stim_duty_cyc = 25;

params_struc2 = setUpStimuli_modular(params_spec2);

%step2.2: Setting up the CS- trial
params_spec3 = params_spec1;
params_spec3.duration = 10;
params_spec3.isi = params_spec3.duration + 120;
if od_to_pair == 3
    params_spec3.odours = 11;
else
    params_spec3.odours = 3;
end
params_spec3.reps = 1;
params_struc3 = setUpStimuli_modular(params_spec3);

%step2.3 Setting up repeats of pairing and CS- trials, n_pairings times
pairing_struc = append_params(params_struc2, params_struc3, 0);
pairing_struc_final = pairing_struc;

for pairing_n = 1:(n_pairings - 1)
    pairing_struc_final = append_params(pairing_struc_final, pairing_struc, 0);
end

params_struc = append_params(params_struc1, pairing_struc_final, 0);


%step3: Adding on the post trials, identical to pre trials
params_struc = append_params(params_struc, params_struc1, 0);
params = params_struc;
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
setup_odor_habituation_trials(curr_dir, 1);
