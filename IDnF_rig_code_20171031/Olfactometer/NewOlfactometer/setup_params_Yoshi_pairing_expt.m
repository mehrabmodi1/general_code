function [params_spec1] = setup_params_Yoshi_pairing_expt(paired_odours, pairing_dur)
%This function sets up a detailed stimulus specification structure and saves
%it into curr_aq_direc to set up stimulus delivery for a pre, pairing and
%post experiment protocol. Paired odor is the vector of odor numbers to be paired.
%if led_elec is 0, led is paired, if led_elec is 1, elec is paired.

%step1:Setting up the pre trials
%creating param specification structure for std odpr resp expt. This is for pre-pairing characterisation. 
params_spec1 = set_ADO_params_default;
params_spec1.odours = [1, 4];
params_spec1.reps = 5; 
params_spec1.duration = 10;
params_spec1.isi = 70;

params_struc1 = setUpStimuli_trains_flex(params_spec1);         %detailed, trial-by-trial parameter specification structure.

params = params_struc1;           %this is the main parameter structure that will ultimately be saved in curr_aq_direc.


%step2: Setting up the pairing trial
%editing param specification structure for only pairing trials. 
params_spec2 = set_ADO_params_LED_stim_only;
params_spec2.reps = 1;
params_spec2.duration = 60;
params_spec2. odours = paired_odours;
params_spec2.led_odours = paired_odours;



%detailed param matrix for a single paired odor trial for CS+
params_paired_tr = setUpStimuli_trains_flex(params_spec2);    

params_spec2unpr = params_spec2;
unpaired_odour = [1, 4];
unpaired_odour(unpaired_odour == paired_odours) = [];
params_spec2unpr.odours = unpaired_odour;
params_spec2unpr.led_odours = [];

%detailed param matrix for a single unpaired odor trial for CS-
params_unpaired_tr = setUpStimuli_trains_flex(params_spec2unpr);

%concatenating pairing trials to params
params = append_params(params, params_paired_tr);
params = append_params(params, params_unpaired_tr);

for p_trial_n = 1:9
    params = append_params(params, params_paired_tr);
    params = append_params(params, params_unpaired_tr);
   
end

%step3: Adding on the post trials
params = append_params(params, params_struc1);
params_mat = params;

%saving params to current acquisition directory
curr_dir = curr_aq_direc;
save([curr_dir, 'params.mat'], 'params_mat');
disp(['saved detailed stim params structure in ', curr_dir]);

