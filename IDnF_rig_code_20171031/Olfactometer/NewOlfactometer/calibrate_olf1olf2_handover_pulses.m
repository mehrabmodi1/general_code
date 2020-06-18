function [params_spec] = calibrate_olf1olf2_handover_pulses(olf1_od, olf2_od)
%This function delivers olf1_od on olf1 and olf2_od on olf2 such that the olf1
%odor pulse is handed over to the olf2 odor pulse. This allows the user to
%observe this handover and adjust the olf1olf2 delay period as necessary to
%make the handover seamless. olf1_od and olf2_od can be vectors of equal
%length and the program will deliver two repeats of pairs of odor pulses
%for each odor pair specified in the two vectors.

params_spec = set_ADO_params_default;

%loop to go through each odor pair
params_spec.reps = 1;
params_spec.odours = olf1_od;
params_spec.odours_olf2 = olf2_od;
params_spec.duration = 10;
params_spec.duration_olf2 = 10;
params_spec.rel_stimLatency_olf2 = 10;
params_spec.isi = 60;
params_spec.post_od_scan_dur = 2;
curr_path = curr_aq_direc;
if exist([curr_path, 'PID_trace_tr-1.mat']) == 2
    delete([curr_path, 'PID_trace_tr-1.mat']);
else
end

present_odours_modular_stim_led(params_spec, 0)

%loading and plotting PID trace

curr_trs = 1;
traces = get_PID_traces(curr_path, curr_trs, 0.1, 1);
figure(1)
plot(traces);
set_xlabels_time(1, 0.1, 10);

%loading in matrix of manually calibrated olf1_olf2 delays for various odor pairs
olf1_olf2_delay_mat = load('E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\olf1_olf2_delay_mat.mat');
olf1_olf2_delay_mat = olf1_olf2_delay_mat.olf1_olf2_delay_mat;
olf1_olf2_delay = olf1_olf2_delay_mat(olf1_od, olf2_od);        %old delay for this odor pair
if isnan(olf1_olf2_delay) == 1
    olf1_olf2_delay = olf1_olf2_delay_mat(10, 3);
else
end

olf1_olf2_delay = input(['Delay currently ', num2str(olf1_olf2_delay), 's. Enter new delay:']);
olf1_olf2_delay_mat(olf1_od, olf2_od) = olf1_olf2_delay;
save('E:\Turner lab\Bitbucket_repos\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\olf1_olf2_delay_mat.mat', 'olf1_olf2_delay_mat');
