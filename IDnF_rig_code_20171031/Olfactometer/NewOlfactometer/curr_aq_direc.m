function [direc, m_direc] = curr_aq_direc()
%This function is for the custom-built setup, it is to be run in the 32-bit
%matlab so as to automatically identify the directory currently being used
%for acquisition. It looks in the main, parent acquisition directory for
%the date-format folder with the most recent (greatest) date number.
%%
mother_direc = 'E:\Turner lab\Data\mehrab\';

% ====  to solve the problem that data are saved in non-mehrab folder ====
%identifying last modified user directory within mother directory
% added by Yichun, Sep 2018
mother_direc = 'E:\Turner lab\Data\';
d = dir(mother_direc);
[~,idx] = max([d.datenum]);
foldername = d(idx).name;
mother_direc = [mother_direc,foldername,'\'];
% ========================================================================

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

%mother directory or directory of the day's expt
m_direc = [mother_direc int2str(max_dir) '\'];


%identifying last modified directory within mother directory
dir_contents = dir(m_direc);
date_time_list = {dir_contents.date};
time_stamp_vec = zeros(1, length(date_time_list));

for item_n = 1:length(date_time_list)
    curr_t = date_time_list{1, item_n};
    time_stamp_vec(1, item_n) = datenum(curr_t(12:20), 'HH:MM:SS');
    
end
time_stamp_vec(1:2) = [];
direc_tail = {dir_contents.name};
direc_tail = {direc_tail{1, 3:length(direc_tail)} };
%finding most recent time stamp
dir_true = 0;
dir_list = {dir_contents.isdir};
dir_list = {dir_list{1, 3:length(dir_list)} };

while dir_true == 0
    [del maxi] = max(time_stamp_vec);
    
    test = dir_list{maxi};
    if test == 1
        dir_true = 1;
    else
        time_stamp_vec(maxi) = 0;
    end
    
end


direc_tail = direc_tail{maxi};

direc = [m_direc direc_tail '\'];

end