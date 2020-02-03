clear all
close all

curr_dir = 'C:\Data\Data\Raw_data\20200131\handover_PID_traces\';

[del, odor_names1] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
[del, odor_names2] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2.xls', 1);
odor_names2{3} = 'Butyl acetate';

[stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig] = load_params_trains_modular(curr_dir, []);    %reading in trial stimulus parameters after matching time stamps to F traces
PA_color = color_vec(2, :);
BA_color = color_vec(1, :);
EL_color = color_vec(3, :);

%identifying stim_mat_simple col numbers
led_on_col_n = find_stim_mat_simple_col('led_on', column_heads);            %identifying relevant column number in stim_mat_simple
od_olf1_col_n = find_stim_mat_simple_col('odor_n', column_heads);           %identifying relevant column number in stim_mat_simple
od_olf2_col_n = find_stim_mat_simple_col('odour_olf2', column_heads);       %identifying relevant column number in stim_mat_simple
dur_olf1_col_n = find_stim_mat_simple_col('duration', column_heads);        %identifying relevant column number in stim_mat_simple
dur_olf2_col_n = find_stim_mat_simple_col('duration_olf2', column_heads);   %identifying relevant column number in stim_mat_simple
od_col_ns = [od_olf1_col_n, od_olf2_col_n];
dur_col_ns = [dur_olf1_col_n, dur_olf2_col_n];
od_durs = unique(stim_mat_simple(:, dur_col_ns(2)));
od_durs(isnan(od_durs)) = [];
odn_list_olf2 = unique(stim_mat_simple(:, od_col_ns(2)));
odn_list_olf2(isnan(odn_list_olf2)) = [];

%plotting simple odor pulse PID traces
%PA
plot_simple_traces(3, 1, PA_color, 1, 2, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir)

%BA
plot_simple_traces(10, 3, BA_color, 3, 4, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir)

%EL
plot_simple_traces(11, 4, EL_color, 5, 6, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir)

%plotting handover odor pulse PID traces
%PA-BA
plot_hover_traces(3, 3, PA_color, BA_color, 7, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir)
%PA-EL
plot_hover_traces(3, 4, PA_color, EL_color, 8, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir)

%BA-PA
plot_hover_traces(10, 1, PA_color, BA_color, 9, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir)
%BA-EL
plot_hover_traces(10, 4, PA_color, EL_color, 10, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir)

%EL-PA
plot_hover_traces(11, 1, EL_color, PA_color, 11, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir)
%EL-BA
plot_hover_traces(11, 3, EL_color, BA_color, 12, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir)




function [] = plot_simple_traces(olf1_odn, olf2_odn, od_color, fign1, fign2, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir)
frame_time = 0.099;

%list of all handover trials, with any combination of odors
handover_trs = find(stim_mat_simple(:, dur_col_ns(1)) == 10 & stim_mat_simple(:, dur_col_ns(2)) == 10);

%plotting olf1 traces
curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == olf1_odn & stim_mat_simple(:, dur_col_ns(1)) == 10);
[del1, del] = intersect(curr_trs, handover_trs);
curr_trs(del) = [];
curr_traces = get_PID_traces(curr_dir, curr_trs, frame_time);
curr_traces = curr_traces(1:(end - 5), 1:2:end);  %getting rid of LED traces
stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
stim_frs = [stim_frs{1}]; 
figure(fign1)
plot(curr_traces, 'lineWidth', 2, 'Color', [0.65, 0.65, 0.65]);
ylabel('PID signal (V)');
set_xlabels_time(fign1, frame_time, 10);
fig_wrapup(fign1, []);
add_stim_bar(fign1, stim_frs, od_color);

%plotting olf2 traces
curr_trs = find(stim_mat_simple(:, od_col_ns(2)) == olf2_odn & stim_mat_simple(:, dur_col_ns(2)) == 10);
[del1, del] = intersect(curr_trs, handover_trs);
curr_trs(del) = [];
curr_traces = get_PID_traces(curr_dir, curr_trs, frame_time);
curr_traces = curr_traces(:, 1:2:end);  %getting rid of LED traces
stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
stim_frs = [stim_frs{2}]; 
figure(fign2)
plot(curr_traces, 'lineWidth', 2, 'Color', [0.65, 0.65, 0.65]);
ylabel('PID signal (V)');
%title ('Pentyl acetate, olfactometer1')
set_xlabels_time(fign2, frame_time, 10);
fig_wrapup(fign2, []);
add_stim_bar(fign2, stim_frs, od_color);

end


function [] = plot_hover_traces(olf1_odn, olf2_odn, od1_color, od2_color, fign, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir)
frame_time = 0.099;

%list of all handover trials, with any combination of odors
handover_trs = find(stim_mat_simple(:, dur_col_ns(1)) == 10 & stim_mat_simple(:, dur_col_ns(2)) == 10);

curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == olf1_odn & stim_mat_simple(:, dur_col_ns(1)) == 10 ...
    & stim_mat_simple(:, od_col_ns(2)) == olf2_odn & stim_mat_simple(:, dur_col_ns(2)) == 10);

curr_traces = get_PID_traces(curr_dir, curr_trs, frame_time);
curr_traces = curr_traces(:, 1:2:end);  %getting rid of LED traces
stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
stim_frs = [stim_frs{1}; stim_frs{2}]; 
figure(fign)
plot(curr_traces, 'lineWidth', 2, 'Color', [0.65, 0.65, 0.65]);
ylabel('PID signal (V)');
%title ('Pentyl acetate, olfactometer1')
set_xlabels_time(fign, frame_time, 10);
fig_wrapup(fign, []);
add_stim_bar(fign, stim_frs, [od1_color; od2_color]);
end


