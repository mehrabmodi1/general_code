function [params_spec1] = setup_params_pairing_expt(paired_odours, led_elec)
%This function sets up a detailed stimulus specification structure and saves
%it into curr_aq_direc to set up stimulus delivery for a pre, pairing and
%post experiment protocol. Paired odor is the vector of odor numbers to be paired.
%if led_elec is 1, led is paired, if led_elec is 2, elec is paired.

%step1:Setting up the pre trials
%creating param specification structure for std KC recording expt. This is for pre-pairing characterisation. 
params_spec1 = set_ADO_params_Yoshi_PA_BA_EL;
params_struc1 = setUpStimuli_trains_flex(params_spec1);         %detailed, trial-by-trial parameter specification structure.
params = params_struc1;           %this is the main parameter structure that will ultimately be saved in curr_aq_direc.


%step2: Setting up the pairing trial
%editing param specification structure for only pairing trials. 
params_spec2 = params_spec1;
params_spec2.duration = 60;
params_spec2.reps = 1;
params_spec2.odours = paired_odours;

if led_elec == 1
    params_spec2.led_odours = paired_odours;
elseif led_elec == 2
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

