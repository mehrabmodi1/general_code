function [path_out] = update_list_path(path_in)
%[curr_dir_list_path] = update_list_path(curr_dir_list_path)
%updating list path to make it suiteable for Brillat

old_i = findstr(path_in, 'code_old');
path_out = [path_in(1:(old_i + 3)), path_in((old_i + 8):end)];
