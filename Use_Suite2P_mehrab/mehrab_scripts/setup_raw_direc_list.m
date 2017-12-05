function [raw_direc_list] = setup_raw_direc_list(raw_direc_base)
prev_direc = pwd;
cd(raw_direc_base);
[~, mat_file_list] = system('dir /S/B params*.*');
mat_file_list = textscan(mat_file_list, '%s', 'Delimiter', '\n');       %list of all param files made by recursively scanning all subdirectories of raw_direc_base
mat_file_list = mat_file_list{1, 1};

n_files = size(mat_file_list, 1);
raw_direc_list = [];
%building raw_direc_list
for file_n = 1:n_files
    curr_direc_str = mat_file_list{file_n, 1};
    dir_starti = findstr(curr_direc_str, '\20');
    dir_endi = findstr(curr_direc_str, '\params_');
    curr_raw_direc = curr_direc_str( (dir_starti + 1):dir_endi);
    
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
    if tot_size > thresh_size
        %making sure that folder structure meets Suite2P requirements (data should be three subfolders away from raw_direc_base)
        slashi = findstr(curr_raw_direc, '\');
        
        while size(slashi, 2) < 3
            if size(slashi, 2) < 3
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
        
        %checking for manual skip_direc tag
        if exist([raw_direc_base, curr_raw_direc, 'skip_direc.txt']) ~=2
            raw_direc_list = [raw_direc_list; {curr_raw_direc}]; 
            
        elseif exist([raw_direc_base, curr_raw_direc, 'skip_direc.txt']) ==2
            
        end
        
        
    
    else
        continue
    end
    
    disp(['building raw data directory list, length is now ' int2str(length(raw_direc_list))])
end


cd(prev_direc);

end