clear all
close all

direc = 'D:\Data\CSHL\20160415\analog_acq\';
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');

data_mat = load([direc, 'acqn_1.mat']);

data_mat = data_mat.data_mat;
dir_contents = dir(direc);
params = [];
for file_n = 1:size(dir_contents, 1)
    curr_name = dir_contents(file_n).name;
    is_param = findstr(curr_name, 'params');
    
    if isempty(is_param) == 0
        if isempty(params) == 1
            params = load([direc curr_name]);
        else
            disp('WARNING: more than one params file found in directory')
        end
    else
    end
    
end

params = params.params;
n_trials = size(params.duration);
odor_t_list = params.odours;
odor_list = unique(odor_t_list);

for odor_n = 1:length(odor_list)
    odor_ni = odor_list(odor_n);
    curr_color = color_vec(odor_ni, :);
    curr_od_trs = find(odor_t_list == odor_ni);
    curr_od_traces = data_mat(:, curr_od_trs);
    
    figure(1)
    plot(curr_od_traces, 'LineWidth', 2, 'Color', curr_color);
    set_xlabels_time(1, 1e-3, 0.5)
    ylabel('PID signal (V)')
    hold on
end

