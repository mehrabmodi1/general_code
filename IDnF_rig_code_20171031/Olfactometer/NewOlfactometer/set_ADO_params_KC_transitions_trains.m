function [params_mat] =  set_ADO_params_KC_transitions_trains(n_reps)

log_file_path = 'E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\reciprocal_expt_logs\KC_transitions_trains.mat';
last_olf1_od = load(log_file_path); 
last_olf1_od = last_olf1_od.last_olf1_od;
last_olf1_od = 3;
odors_to_switch = [3, 11];     %PA and EL on olf1
od_remi = find(odors_to_switch == intersect(last_olf1_od, odors_to_switch));
odors_to_switch(od_remi) = [];
olf1_od = odors_to_switch;
if olf1_od == 3 %PA on olf1
    olf2_od  = 4; %EL on olf2
elseif olf1_od == 11 %EL on olf1
    olf2_od = 1; %PA on olf2
else
end

last_olf1_od = olf1_od;
save(log_file_path, 'last_olf1_od');


%setting up simple trials for olf1
params_spec1 = set_ADO_params_default;
params_spec1.duration = [10, 30];
params_spec1.reps = n_reps;
params_spec1.post_od_scan_dur = 20;
params_spec1.isi = 80;

params_spec1.odours = olf1_od;
params_spec1.odours_olf2 = [];
params_spec1.duration_olf2 = [10, 30];
params_struc1 = setUpStimuli_modular(params_spec1); 

params_spec1_2 = params_spec1;
params_spec1_2.duration = [0.1, 0.11];
params_spec1_2.odours_olf2 = olf2_od;
params_spec1_2.duration_olf2 = [10, 30];
params_struc1_2 = setUpStimuli_modular(params_spec1_2); 

params_struc = append_params(params_struc1, params_struc1_2, 1);  %combining and randomising explicit param spec structures

%setting up handover trials olf1_od to olf2_od
params_spec2 = params_spec1;
params_spec2.reps = 3;
params_spec2.odours_olf2 = olf2_od;
params_spec2.duration = 10;
params_spec2.duration_olf2 = 10;
params_spec2.rel_stimLatency_olf2 = 10;
params_struc2 = setUpStimuli_modular(params_spec2);
%params_struc = append_params(params_struc, params_struc2, 1);  %combining and randomising explicit param spec structures


%setting up handover trials olf2_od to olf1_od
params_spec2_2 = params_spec2;
if olf1_od == 3 %PA on olf1
    params_spec2_2.odours_olf2 = 1; %PA on olf2
    params_spec2_2.odours = 11; %EA on olf
elseif olf1_od == 11 %EA on olf1
    params_spec2_2.odours_olf2 = 4; %EA on olf2
    params_spec2_2.odours = 3; %PA on olf
else
end
params_struc2_2 = setUpStimuli_modular(params_spec2_2);
%params_struc = append_params(params_struc, params_struc2_2, 1);  %combining and randomising explicit param spec structures


%setting up rand train trials olf1
params_spec3 = params_spec1;
params_spec3.reps = 3;
params_spec3.isi = 120;
params_spec3.duration = 60;
params_spec3.rand_trains = 1;
params_spec3.n_rand_trains = 1;
params_spec3.mean_rand_pulse_dur = 6;
params_struc3 = setUpStimuli_modular(params_spec3);
%getting rid of embedded simple trials
for tr_n = size(params_struc3, 2):-1:1      %going in reverse to prevent frame-shifts below
    if size(params_struc3(tr_n).pulse_train, 1) == 1
        params_struc3(tr_n) = [];
    else
    end
end
%params_struc = append_params(params_struc, params_struc3, 1);  %combining and randomising explicit param spec structures
train_olf1 = params_struc3(1).pulse_train;


%setting up rand train trials olf2
params_spec3_2 = params_spec1_2;
params_spec3_2.reps = n_reps;
params_spec3_2.duration = [0.1];
params_spec3_2.duration_olf2 = 60;
params_spec3_2.rand_trains_olf2 = 1;
params_spec3_2.n_rand_trains_olf2 = 1;
params_spec3_2.mean_rand_pulse_dur_olf2 = 1;
params_struc3_2 = setUpStimuli_modular(params_spec3_2);
%getting rid of embedded simple trials
for tr_n = size(params_struc3_2, 2):-1:1    %going in reverse to prevent frame-shifts below
    if size(params_struc3_2(tr_n).pulse_train_olf2, 1) == 1
        params_struc3_2(tr_n) = [];
    else
    end
end
%params_struc = append_params(params_struc, params_struc3_2, 1);  %combining and randomising explicit param spec structures

%setting up combined, olf1, olf2 rand_train trials
params_struc4 = params_struc3_2;
for tr_n = 1:params_spec3_2.reps
    params_struc4(tr_n).duration = 60;
    params_struc4(tr_n).pulse_train = train_olf1;
end
params_struc = append_params(params_struc, params_struc4, 1);  %combining and randomising explicit param spec structures
params_mat = params_struc;

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

setup_odor_habituation_trials(curr_dir, 0);

