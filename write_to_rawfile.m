function [] = write_to_rawfile(dataset, trial_n, new_matrix)
%This function over-writes the raw data matrix associated with a 2P object
%(dataset) for a given trial (trial_n) with the data in 'new_matrix'.

% dataset = load('C:\Data\CSHL\test_movt - Copy\expt\expt.mat');
% dataset = dataset.data;
% trial_n = 1;
% new_matrix = load('C:\Data\CSHL\test_movt - Copy\expt\test_movt\rawData1.mat');
% new_matrix = new_matrix.rawData1;
% new_matrix(80:85, 80:85, 1:100) = 0.25;
% keyboard
eval(['rawData' int2str(trial_n) '= new_matrix;' ])

%folder within which raw data file is to be re-written
write_direc = dataset(trial_n).info.rawDataDir;

%over-writing raw data file with new_matrix
save([write_direc '\rawData' int2str(trial_n) '.mat'], ['rawData' int2str(trial_n)] );



end