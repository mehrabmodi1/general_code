function [direc] = curr_aq_direc()
%This function is for the custom-built setup, it is to be run in the 32-bit
%matlab so as to automatically identify the directory currently being used
%for acquisition. It looks in the main, parent acquisition directory for
%the date-format folder with the most recent (greatest) date number.


mother_direc = 'C:\Data\CSHL\';         %IMPORTANT - set this to the parent directory within which to look for current acquisition directory!

dir_contents = dir(mother_direc);
dir_contents = {dir_contents.name};

max_dir = 0;
for dir_count = 1:length(dir_contents)
    curr_dir = dir_contents{1, dir_count};
    if isdir([mother_direc curr_dir]) == 0
        continue
    else
    end
    curr_dir = str2num(curr_dir);
    if isempty(curr_dir) == 1
        continue
    else
    end
    max_dir = max([max_dir, curr_dir]);
    
    
end

direc = [mother_direc int2str(max_dir) '\'];

end