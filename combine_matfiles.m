function [data] = combine_matfiles(direc)
%This function reads in the separate mat files generated for the datasets
%from a given experiment folder. The combined mat file is saved in the expt
%folder.

if nargin < 1
    direc = pwd;
else
end

dir_contents = dir('*.mat');
d = dir_contents;
d(find(strcmp({d.name},'.')))=[];
d(find(strcmp({d.name},'..')))=[];
d(find(strcmp({d.name},'expt.mat')))=[];
d(find(strcmp({d.name},'expt_raw_traces.mat')))=[];
d(find(strcmp({d.name},'aved_tr_frames.mat')))=[];
d(find(strcmp({d.name},'movie.mat')))=[]
dir_contents = d;

n_sets = size(dir_contents, 1);

%all_mats = {n_sets, 1};
%timestamps = zeros(n_sets, 1);
set_counter = 0;

for set_n = 1:n_sets
    
    curr_fname = dir_contents(set_n).name;
    %skipping .mat files that aren't part of this procedure
    if strcmp(curr_fname, 'expt.mat') == 1
        continue
    elseif strcmp(curr_fname, 'expt_raw_traces.mat') == 1
        continue
    else
    end
    set_counter = set_counter + 1;
    curr_mat = load([direc '\' curr_fname]);
    
    curr_mat = curr_mat.data;
    all_mats{set_counter, 1} = curr_mat;
    try
        timestamps(set_counter, 1) = curr_mat(1).stim.timestamp;
    catch
        if n_sets == 1
            timestamps(set_counter, 1) = 1;
        else
            beep
            keyboard
        end
    end
    
end
n_sets = set_counter;
clear set_counter;

%sorting set numbers by time of acquisition before combining them into a single mat file
timestamps = [timestamps, (1:n_sets)'];
timestamps = sortrows(timestamps);
n_trials_vec = [];

for set_n = 1:n_sets
    set_ns = timestamps(set_n, 2);
    curr_mat = all_mats{set_ns};
    
    if set_n == 1
        data = curr_mat; 
    else
    end
    
    if set_n > 1
        n_start = size(data, 2) + 1;
        n_end = size(data, 2) + size(curr_mat, 2);
        try
            data(n_start:n_end) = curr_mat(1:size(curr_mat, 2));
        catch
            keyboard
        end
    else
    end
    n_trials = size(curr_mat, 2);
    n_trials_vec = [n_trials_vec; n_trials];
    
end

save([direc '\trial_n_list.txt'], 'n_trials_vec', '-ASCII');

end