function [params_spec] = calibrate_olf1olf2_handover_pulses(olf1_od, olf2_od)
%This function delivers olf1_od on olf1 and olf2_od on olf2 such that the olf1
%odor pulse is handed over to the olf2 odor pulse. This allows the user to
%observe this handover and adjust the olf1olf2 delay period as necessary to
%make the handover seamless. olf1_od and olf2_od can be vectors of equal
%length and the program will deliver two repeats of pairs of odor pulses
%for each odor pair specified in the two vectors.

params_spec = set_ADO_params_default;

%loop to go through each odor pair
params_spec.reps = 2;
params_spec.odours = olf1_od;
params_spec.odours_olf2 = olf2_od;
params_spec.duration = 10;
params_spec.duration_olf2 = 10;
params_spec.rel_stimLatency_olf2 = 10;
params_spec.isi = 60;

