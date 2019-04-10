function [dataset_ch] = change_tpobj_direc(direc) 
    %This function edits the folder path in an existing 2P object in case it
    %has been transferred from one path to another. It then re-saves the 2P
    %object.

%testing lines
%     direc = 'D:\Data\CSHL\20150917\fly2\';
    

    dataset = load([direc 'expt.mat']);
    dataset = dataset.data;


    %fixing path in two photon object file
    obj_dir = dataset(1).info.Filename;
    if obj_dir(1) == 'C'    %test for old, incorrect path
        curr_dir_r = dataset(1).info.rawDataDir;
        curr_dir_r(1) = 'D';

        %looping through each trial in the dataset
        for tr_n = 1:size(dataset, 2)
            curr_dir = dataset(tr_n).info.Filename;
            curr_dir(1) = 'D';

            dataset(tr_n).info.rawDataDir = curr_dir_r;
            dataset(tr_n).info.Filename = curr_dir;

        end
        data = dataset;
        save([direc 'expt.mat'], 'data');
        clear data
    else
    end

end