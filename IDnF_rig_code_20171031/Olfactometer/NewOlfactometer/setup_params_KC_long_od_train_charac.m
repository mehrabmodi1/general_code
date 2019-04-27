function [params_spec1] = setup_params_KC_long_od_train_charac()
%syntax: function [params_spec1] = setup_params_pairing_expt(paired_odours, led_elec, CS_dur, US_dur, US_dutycyc)
%This function sets up a detailed stimulus specification structure and saves
%it into curr_aq_direc to set up stimulus delivery for long odor followed by rapid train protocols. 

%step1:Setting up the long od trials
%creating param specification structure for std KC recording expt. This is for pre-pairing characterisation. 
params_spec1 = set_ADO_params_charac_KC_long_resps;
params_struc1 = setUpStimuli_trains_flex(params_spec1);         %detailed, trial-by-trial parameter specification structure.
params = params_struc1;           %this is the main parameter structure that will ultimately be saved in curr_aq_direc.


%step2: Setting up the foreground train trials
%editing param specification structure for trains. 
params_spec2 = params_spec1;
params_spec2.duration = 6;
params_spec2.n_od_pulses = 10;
params_spec2.inter_pulse_interval = 21;
params_spec2.reps = 3;
params_spec2.isi = 250;
params_spec2.pulse_type = 1;

params_struc2 = setUpStimuli_trains_flex(params_spec2);         %detailed, trial-by-trial parameter specification structure.
%concatenating train trials to params
params = append_params(params, params_struc2);


%step3: Setting up the background train trials
%editing param specification structure for trains. 
params_spec3 = params_spec1;
params_spec3.duration = 2;
params_spec3.n_od_pulses = 10;
params_spec3.inter_pulse_interval = 21;
params_spec3.reps = 3;
params_spec3.isi = 250;
params_spec3.pulse_type = 1;

params_struc3 = setUpStimuli_trains_flex(params_spec3);         %detailed, trial-by-trial parameter specification structure.
%concatenating train trials to params
params = append_params(params, params_struc3);
params_mat = params;

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

