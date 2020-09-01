function [params_spec1] =  set_ADO_params_alpha1_post_ar_training()


%setting up simple trials for olf1
params_spec1 = set_ADO_params_default;

params_spec1.duration = [10];
params_spec1.reps = 7;
params_spec1.post_od_scan_dur = 20;
params_spec1.isi = 60;
params_spec1.odours = [3, 9, 11];
params_spec1.odours_olf2 = [];