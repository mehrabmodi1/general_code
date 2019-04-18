function [data_mat] = load_suite2P_results(direc)
%This function identifies the most recently created Suite2P results file
%(such as the one created after manually pruning ROIs) and loads it into
%memory.

prev_direc = pwd;
cd(direc);
dir_contents = dir;
dir_contents(1:2) = [];
n_files = size(dir_contents, 1);
%finding most recent file
max_datenum = [];
for file_n = 1:n_files
    fname = dir_contents(file_n).name;
    if isempty(findstr(fname, 'Nk200')) == 1
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
disp(['Loading Suite2P results file ' dir_contents(max_datenum(2)).name])
data_mat = load([direc, dir_contents(max_datenum(2)).name]);
try
    data_mat = data_mat.dat;
catch
end
disp('Done Loading.')
cd(prev_direc);
end