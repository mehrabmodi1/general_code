function [] = tune_olf1(od_n)

params = set_ADO_params_default;

params.odours = od_n;
params.reps = 1;
params.isi = 2500;
params.duration = 10;
params.n_od_pulses = 50;
params.inter_pulse_interval = 10;
params.odours_olf2 = [];

present_odours_modular_stim_led(params, 0)
