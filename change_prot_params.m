clear all
close all

list_direc = ['D:\Data\CSHL\dataset_list_glc_ara_20150917.txt'];
fid = fopen(list_direc);

direc_counter = 0;

%loop to go through all experiment datasets listed in list file
while 1
    direc_counter = direc_counter + 1;
    direc = fgetl(fid);
    %replacing C: with D:
    a = findstr(direc, 'C:');
    direc(a) = 'D';
    
       
    if ischar(direc) ~= 1
        break
    else
    end
    
    direc = [direc '\'];
    
    cd(direc)
    
    dataset = load([direc 'expt.mat']);
    dataset = dataset.data;
    
    n_trials = size(dataset, 2);
    
    for trial_n = 1:n_trials
        dataset(trial_n).stim.prot_switch_trials = ['18'; '36'];
    end
    
    curr_set = 'expt.mat';
    
    data = dataset;
    save([direc curr_set], 'data');
   
end
fclose(fid);