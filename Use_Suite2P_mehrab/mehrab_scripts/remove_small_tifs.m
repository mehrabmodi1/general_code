function [] = remove_small_tifs(raw_direc)
%syntax: [] = remove_small_tifs(raw_direc)
%function to get rid of empty .tif files created by ScanImage sometimes

old_direc = pwd;
cd(raw_direc);
dir_contents = dir('*.tif');

small_files = [];
%detecting small files
for f_num = 1:size(dir_contents, 1)
    if dir_contents(f_num).bytes < 15500000
        small_files = [small_files; f_num];
    else
    end
end
%deleting small files
if isempty(small_files) == 0
    for s_f_num = 1:length(small_files)
        s_f_num_i = small_files(s_f_num);
        curr_fname = dir_contents(s_f_num_i).name;
        
        if isempty(findstr(curr_fname, '.tif')) == 0
            delete(curr_fname);
        else
        end
    end
else
end

cd(old_direc)
end