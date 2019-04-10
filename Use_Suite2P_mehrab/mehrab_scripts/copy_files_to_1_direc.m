%checking if directory structure was extended with a \1\ folder at the
%end and copying over the results file from ROI_prune to that folder

function [] = copy_files_to_1_direc(results_direc, raw_direc)
one_i = findstr(raw_direc, '\1\');
if isempty(one_i) == 0

    raw_direc_t = raw_direc(1:one_i);
    cd([results_direc, raw_direc_t]);
    dir_contents_t = dir('*.mat');

    for file_n_t = 1:size(dir_contents_t, 1)
        f_name_t = dir_contents_t(file_n_t).name;
        copyfile([results_direc, raw_direc_t, f_name_t], [results_direc, raw_direc]);
        delete([results_direc, raw_direc_t, f_name_t]);
    end
    cd([results_direc, raw_direc])
else
end