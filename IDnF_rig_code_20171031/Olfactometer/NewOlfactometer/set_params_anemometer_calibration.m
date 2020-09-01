clear all
close all
%Note: empty odor vials should be connected to odor port 3 on olfactometer1
%and odor port 1 on olfactometer 2.
params = set_ADO_params_default;

%olf1 trials
params.reps = 5;
params.odours = 3;
params.duration = 10;
params.isi = 60;
params.odours_olf2 = [];
params_struc1 = setUpStimuli_modular(params); 


%olf2 trials
params2 = params;
params2.duration = 0.1;
params2.odours_olf2 = 1;
params2.duration_olf2 = 10;
params2.rel_stimLatency_olf2 = 0;

params_struc2 = setUpStimuli_modular(params2); 
params_struc = append_params(params_struc1, params_struc2, 0);

%olf1-olf2 transition trials
params3 = params2;
params3.duration = 10;
params3.rel_stimLatency_olf2 = 10;

params_struc3 = setUpStimuli_modular(params3); 
params_struc = append_params(params_struc, params_struc3, 0);

params_mat = params_struc;
curr_dir = curr_aq_direc;
save([curr_dir, 'params.mat'], 'params_mat');
