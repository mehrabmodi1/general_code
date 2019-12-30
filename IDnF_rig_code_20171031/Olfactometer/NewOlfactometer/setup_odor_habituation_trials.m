function [] = setup_odor_habituation_trials(path)
%This function looks in 'path' for and loads in an explicit stimulus specification
%structure. It then lists the odors being delivered on each olfactometer,
%and creates two, 20s pre-exposure trials for each odor, randomly
%interleaved. Scanning will not be triggered during these trials. They are
%added on at the beginning of the stimulus structure and written to file in
%the same location.

%loading explicit param structure
params_struc_orig = load([path, 'params.mat']);
params_struc_orig = params_struc_orig.params_mat;

if params_struc_orig(1).trigger_scan == 0
    redo = input('Trials with no scan trigger detected, redo habituation trials?, 1 - Yes, 0 - No');
    
    if redo == 0
        return
    else
    end

else
end
    
%building lists of odors delivered on olf1 and olf2
n_trs = size(params_struc_orig, 2);
olf1_ods = zeros(n_trs, 1) + nan;
olf2_ods = zeros(n_trs, 1) + nan;
for tr_n = 1:n_trs
    olf1_ods(tr_n, 1) = params_struc_orig(tr_n).odours;
    olf2_ods(tr_n, 1) = params_struc_orig(tr_n).odours_olf2;
    
end

%getting rid of nans
olf1_ods(isnan(olf1_ods)) = [];
olf2_ods(isnan(olf2_ods)) = [];

olf1_ods = unique(olf1_ods);
olf2_ods = unique(olf2_ods);


%building implicit stim param structure for habituation trials
%for olf1 odors
params_spec1 = set_ADO_params_default;
params_spec1.trigger_scan = 0;
params_spec1.duration = 20;
params_spec1.reps = 2;
params_spec1.odours = olf1_ods;
params_spec1.odours_olf2 = [];
params_spec1.isi = 50;
params_spec1.post_od_scan_dur = 0;
params_spec1.stimLatency = 0;

params_spec1.odourNames = params_struc_orig(1).odourNames;
params_spec1.odourNames_olf2 = params_struc_orig(1).odourNames_olf2;


params_struc1 = setUpStimuli_modular(params_spec1);         %explicit, trial-by-trial parameter specification structure for olf1 odours.
params_struc = params_struc1;

%for olf2 odors
if isempty(olf2_ods) == 0
    params_spec2 = params_spec1;
    params_spec2.odours_olf2 = olf2_ods;
    params_spec2.odours = olf1_ods(1:length(olf2_ods));
    params_spec2.duration = 0.1;            %no odor to be delivered on olf1
    params_spec2.duration_olf2 = 20;
    params_spec2.rel_stimLatency = 0;
    
    params_struc2 = setUpStimuli_modular(params_spec2);         %explicit, trial-by-trial parameter specification structure for olf1 odours.
    params_struc = append_params(params_struc, params_struc2, 1);  %combining and randomising explicit param spec structures

else
end

%adding these habituation trials to the beginning of loaded param
%specification structure
params_mat = append_params(params_struc, params_struc_orig, 0);

%saving to file
save([path, 'params.mat'], 'params_mat');
