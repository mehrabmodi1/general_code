function [params_spec1] = setup_params_pairing_expt_integrator_HepIAA(US_dutycyc, LED_power)
%syntax: function [params_spec1] = setup_params_pairing_expt(paired_odours, led_elec, CS_dur, US_dur, US_dutycyc)
%This function sets up a detailed stimulus specification structure and saves
%it into curr_aq_direc to set up stimulus delivery for a pre, pairing and
%post experiment protocol. Paired odor is the vector of odor numbers to be paired.
%if led_elec is 0, led is paired, if led_elec is 1, elec is paired.

%step1:Setting up the pre trials
%creating param specification structure for std KC recording expt. This is for pre-pairing characterisation. 
params_spec1 = set_ADO_params_HepIAA;
params_spec1.reps = 3;
params_spec1.duration = 10;
params_spec1.isi = 60;
params_spec1.duration_olf2 = 0;
params_struc1 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.

params_spec1.duration_olf2 = 10;
params_spec1.duration = 0.05;
params_struc2 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.
params_struc_pre = append_params(params_struc1, params_struc2, 1);
params = params_struc_pre;           %this is the main parameter structure that will ultimately be saved in curr_aq_direc.



%step2: Setting up the pairing trial
%editing param specification structure for only pairing trial. 
params_spec2 = set_ADO_params_fore_distracter_HepIAA;
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
    elseif ovwrite == 0
    end
else
    save([curr_dir, 'params.mat'], 'params_mat');
end

disp(['saved detailed stim params structure in ', curr_dir]);

