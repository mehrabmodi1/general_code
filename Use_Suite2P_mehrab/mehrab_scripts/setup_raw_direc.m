function [raw_direc] = setup_raw_direc(raw_direc_base, raw_direc)
prev_direc = pwd;
curr_raw_direc = raw_direc;

%checking if setup_raw_direc has already been run on current folder
if isfolder([raw_direc_base, curr_raw_direc, '1\']) == 1
    raw_direc = [curr_raw_direc, '1\'];
    return
else
end

%checking if there are enough frames in current dataset to run Suite2P
cd([raw_direc_base, curr_raw_direc]);
dir_contents = dir('*.tif');
n_tifs = size(dir_contents, 1);
tot_size = 0;                       %keeping track of the total bytes of tif frames as a proxy for number of acquired frames in dataset (must be > 500 to run Suite2P)
for tif_n = 1:n_tifs
    curr_size = dir_contents(tif_n).bytes;
    tot_size = tot_size + curr_size;
end

thresh_size = 500e6;                %at least 500 MB of frames needed before running Suite2P on dataset
if tot_size > thresh_size && n_tifs >= 8
    %making sure that folder structure meets Suite2P requirements (data should be three subfolders away from raw_direc_base)
    slashi = findstr(curr_raw_direc, '\');
    if size(slashi, 2)  < 3
        while size(slashi, 2) < 3
            if size(slashi, 2) < 3
                
                %PICK UP THREAD HERE
                %put in condition to check if \1 folder has already been
                %created and data copied over.
                
                dir_contents_data = dir([raw_direc_base, curr_raw_direc]);
                dir_contents_data(1:2) = [];
                mkdir([raw_direc_base, curr_raw_direc, '\1']);

                n_raw_files = size(dir_contents_data, 1);
                for raw_file_n = 1:n_raw_files
                    curr_raw_name = dir_contents_data(raw_file_n).name;
                    copyfile([raw_direc_base, curr_raw_direc, curr_raw_name], [raw_direc_base, curr_raw_direc, '\1\', curr_raw_name]);
                    delete([raw_direc_base, curr_raw_direc, curr_raw_name]);
                end
                curr_raw_direc = [curr_raw_direc, '1\'];
                slashi = findstr(curr_raw_direc, '\');
                
            else
            end
        end
    else
    end
    raw_direc = curr_raw_direc;
    %checking for manual skip_direc tag
    if exist([raw_direc_base, curr_raw_direc, 'skip_direc.txt']) ==2
        raw_direc = [];         %skipping this direc
    else

    end



else
    raw_direc = [];         %skipping this direc
end
   
cd(prev_direc);

end