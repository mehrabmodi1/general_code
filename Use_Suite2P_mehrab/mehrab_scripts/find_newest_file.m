function [newest_filename] = find_newest_file(direc, sub_string)
%suntax: [newest_filename] = find_newest_file(direc, sub_string)
%This function finds the newest file containing sub_string in its filename
%in the directory direc. Use '.' as sub_string if sub-string specification is 
%not needed.

prev_direc = pwd;
cd(direc)

dir_contents = dir;
dir_contents(1:2) = [];
n_files = size(dir_contents, 1);
%finding most recent file
max_datenum = [];
for file_n = 1:n_files
    fname = dir_contents(file_n).name;
    if isempty(findstr(fname, sub_string)) == 1
        continue
    else
    end
    curr_datenum = dir_contents(file_n).datenum;
    if isempty(max_datenum) == 1
        max_datenum = [curr_datenum, file_n];
    elseif isempty(max_datenum) == 0
        if max_datenum(1) < curr_datenum == 1
            max_datenum(1) = curr_datenum;
            max_datenum(2) = file_n;
        else
        end
    end

end
try
    newest_filename = dir_contents(max_datenum(2)).name;
catch
    keyboard
end
cd(prev_direc)
