function [params_spec1] = setup_params_pairing_expt_bkwd(paired_odours, led_elec, CS_dur, US_dur, US_dutycyc)
%This function sets up a detailed stimulus specification structure and saves
%it into curr_aq_direc to set up stimulus delivery for a pre, pairing and
%post experiment protocol. Paired odor is the vector of odor numbers to be paired.
%if led_elec is 0, led is paired, if led_elec is 1, elec is paired.

%step1:Setting up the pre trials
%creating param specification structure for std KC recording expt. This is for pre-pairing characterisation. 
params_spec1 = set_ADO_params_Yoshi_PA_BA_EL;
params_struc1 = setUpStimuli_trains_flex(params_spec1);         %detailed, trial-by-trial parameter specification structure.
params = params_struc1;           %this is the main parameter structure that will ultimately be saved in curr_aq_direc.


%step2: Setting up the pairing trial
%editing param specification structure for only pairing trials. 
params_spec2 = params_spec1;
params_spec2.duration = CS_dur;
stimLatency_old = params_spec2.stimLatency;
params_spec2.stimLatency = stimLatency_old + US_dur;
params_spec2.isi = params_spec2.isi + US_dur;
params_spec2.reps = 1;
params_spec2.odours = paired_odours;
params_spec2.stim_init_delay_ms = (stimLatency_old.*1000) + 750;     %750 ms added on for odor to flow through tube and reach fly
params_spec2.stim_dur = 1000.*US_dur;                                    %in ms
params_spec2.st_duty_cyc = US_dutycyc;                                                 % reduced LED stim duty cycle to 5% form 50%

if led_elec == 0
    params_spec2.led_odours = paired_odours;
elseif led_elec == 1
    params_spec2.elec_odours = paired_odours;
else
end
params_struc2 = setUpStimuli_trains_flex(params_spec2);         %detailed, trial-by-trial parameter specification structure.
%concatenating pairing trials to params
params = append_params(params, params_struc2);

%adding on trials for CS minus presentation
other_ods = params_spec1.odours;
other_ods(other_ods == paired_odours) = [];
params_spec2.odours = other_ods;
params_spec2.led_odours = [];
params_spec2.elec_odours = [];
params_struc2 = setUpStimuli_trains_flex(params_spec2);         %detailed, trial-by-trial parameter specification structure.
%concatenating CS minus trials to params
params = append_params(params, params_struc2);


%step3: Adding on the post trials
params = append_params(params, params_struc1);
params_mat = params;

%saving params to current acquisition directory
curr_dir = curr_aq_direc;
save([curr_dir, 'params.mat'], 'params_mat');
disp(['saved detailed stim params structure in ', curr_dir]);

