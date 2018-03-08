function [dir_contents_sorted] = dir_date_sorted(direc, dir_string)
%This function reads in the contents of the folder direc, using the syntax
%dir('dir_string') after making direc the current folder (see documentation 
%for dir() for how this could be useful). If you don't want to specify a 
%search string, set dir_string = []. The function then sorts the directory 
%contents by their datenums.
%Mehrab Modi 20180304

%Test variables
% direc = 'C:\Data\Code\general_code\test';
% dir_string = [];

old_direc = pwd;
cd(direc);

if isempty(dir_string) == 0
    dir_contents = dir(dir_string);
else
    dir_contents = dir();
    dir_contents(1:2) = [];
end
n_files = size(dir_contents, 1);
date_list = zeros(n_files, 2) + nan;

%building list of datenums
for file_n = 1:n_files          %Note: This will also sort folders, if not excluded by the dir search string
    date_list(file_n, 1) = dir_contents(file_n).datenum;
end

%sorting list of datenums
date_list(:, 2) = (1:1:n_files)';
date_list = sortrows(date_list);

%sorting list of directory contents
for file_n = 1:n_files
    orig_list_n = date_list(file_n, 2);
    dir_contents_sorted(file_n, 1) = dir_contents(orig_list_n, 1);
end

cd(old_direc)