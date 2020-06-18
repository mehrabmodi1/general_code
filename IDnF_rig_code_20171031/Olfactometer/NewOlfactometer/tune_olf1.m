function [] = tune_olf1(od_n, mode)

params = set_ADO_params_default;

params.odours = od_n;
params.odours_olf2 = [];

if strcmp(mode, 'n') == 1   %rapidly deliver many pulses to tune olfactometer
   
    params.reps = 1;
    params.isi = 2500;
    params.duration = 10;
    params.n_od_pulses = 50;
    params.inter_pulse_interval = 10;
    
    %delivering many pulses in quick succession
    present_odours_modular_stim_led(params, 0)
    
elseif strcmp(mode, 'plot') == 1
    params.reps = 1;
    params.isi = 35;
    params.duration = 10;
    params.n_od_pulses = 1;
    params.inter_pulse_interval = 10;
    params.post_od_scan_dur = 2;
    
    %delivering pulses as trials to aquire multiple PID traces
    present_odours_modular_stim_led(params, 0)
    
    %loading and plotting PID traces
    curr_path = curr_aq_direc;
    curr_trs = 1:1:5;
    traces = get_PID_traces(curr_path, curr_trs, 0.1, 1);
    
    plot(traces)
    
else
end

