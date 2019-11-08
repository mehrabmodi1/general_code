function [params] = set_params_single_ax_plasticity(paired_odor, odor_list)

%step1:Setting up the pre trials
%step1.1 setting up pre-post odours to be delivered with olf1
%creating param specification structure for std KC recording expt. This is for pre-pairing characterisation. 
params_spec1 = set_ADO_params_default;

params_spec1.odours = odor_list;
params_spec1.reps =5;
params_spec1.duration = 60;
params_spec1.isi = 120;
params_spec1.duration_olf2 = 0;
params_struc1 = setUpStimuli_modular(params_spec1);         %detailed, trial-by-trial parameter specification structure for olf1 odour.
params = params_struc1;

%step2: Setting up the pairing trial
%editing param specification structure for only pairing trial. 
params_spec2 = params_spec1;
params_spec2.led_odours = paired_odor;
params_spec2.odours = paired_odor;
params_spec2.reps = 1;
params_spec2.duration = 60;
params_spec2.isi = 130;
params_spec2.rel_stim_init_delay = 5;
params_spec2.stim_dur = 60;
params_spec2.stim_freq = 1;
params_spec2.st_duty_cyc = 0.5; %percent


params_struc_pairing = setUpStimuli_modular(params_spec2);         %detailed, trial-by-trial parameter specification structure.

%concatenating pairing trial to params
params = append_params(params, params_struc_pairing, 0);


%step3: Adding on the post trials
params = append_params(params, params_struc1, 0);
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
