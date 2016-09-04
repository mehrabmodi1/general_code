function [isi, odor_t_list, odorNames, n_reps, stim_time] = load_stim_param(direc)
%function to read in stimulus parameters for an odor-evoked activity
%dataset
%Mehrab Modi, 20140926

dir_contents = dir([direc '*.mat']);
fname = dir_contents.name;

params = load([direc fname]);
params = params.params;

isi = params.isi;
odor_t_list = params.odours;
odorNames = params.odourNames;
n_reps = params.reps;
stim_time = params.stimLatency;
end