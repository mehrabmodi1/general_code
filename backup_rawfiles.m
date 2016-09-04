function [] = backup_rawfiles(dataset)
%This function creates copies of the raw data files after rigid-body motion
%correction so that quick reversion is possible.

% dataset = load('C:\Data\CSHL\test_movt - Copy\expt\expt.mat');
% dataset = dataset.data;

backup_direc = dataset(1).info.rawDataDir;
write_direc = [backup_direc '\backup_rawfiles\'];
mkdir(write_direc);

n_trials = size(dataset, 2);

for trial_n = 1:n_trials
    if exist([write_direc '\rawData' int2str(trial_n) '.mat']) == 2
        continue
    else
        curr_matrix = load([backup_direc '\rawData' int2str(trial_n) '.mat']);      %loading rawfile data
        eval(['curr_matrix = curr_matrix.rawData' int2str(trial_n) ';']);

        var_name = genvarname(['rawData' int2str(trial_n)]);
        eval([var_name '= curr_matrix;' ]);                        %generating variable with appropriate name to hold rawfile data

        save([write_direc '\rawData' int2str(trial_n) '.mat'], var_name);
        disp(['backed up raw file ' int2str(trial_n)])
    end
end


end