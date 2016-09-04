function [] = recover_rawfiles(dataset)
% dataset = load('C:\Data\CSHL\test_movt - Copy\expt\expt.mat');
% dataset = dataset.data;

write_direc = dataset(1).info.rawDataDir;
backup_direc = [write_direc '\backup_rawfiles\'];

n_trials = size(dataset, 2);

for trial_n = 1:n_trials
    curr_matrix = load([backup_direc '\rawData' int2str(trial_n) '.mat']);      %loading rawfile data
    eval(['curr_matrix = curr_matrix.rawData' int2str(trial_n) ';']);
    
    var_name = genvarname(['rawData' int2str(trial_n)]);                        %generating variable with appropriate name to hold rawfile data
    eval([var_name '= curr_matrix;' ]);                        
    
    save([write_direc '\rawData' int2str(trial_n) '.mat'], var_name);           %over-writing existing rawfile with backup file
    disp(['recovered raw file ' int2str(trial_n)])
end

end