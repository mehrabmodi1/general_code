clear all
close all


MPPC_trace_path = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\Janelia\Hardware_stuff\SiPM_testing\SiPM_raw_traces.xls';
MPPC_traces = xlsread(MPPC_trace_path, 1);
bk = mean(MPPC_traces(:, 4));
MPPC_traces = MPPC_traces - bk;

bk_vals = mean(MPPC_traces(1:500, :));
bk_mat = repmat(bk_vals, size(MPPC_traces, 1), 1);
dff_MPPC_traces = (MPPC_traces - bk_mat)./bk_mat;


PMT_trace_path = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\Janelia\Hardware_stuff\SiPM_testing\PMT_raw_traces.xls';
PMT_traces = xlsread(PMT_trace_path, 1);
bk = mean(PMT_traces(:, 4));
PMT_traces = PMT_traces - bk;

bk_vals = mean(PMT_traces(1:500, :));
bk_mat = repmat(bk_vals, size(PMT_traces, 1), 1);
dff_PMT_traces = (PMT_traces - bk_mat)./bk_mat;

stim_frs = [774, 1107];

figure(1)
plot(MPPC_traces)
set_xlabels_time(1, 0.03, 10);
fig_wrapup(1)
add_stim_bar(1, stim_frs, [0, 0, 0])

figure(2)
plot(dff_MPPC_traces(:, 1:3))
set_xlabels_time(2, 0.03, 10);
fig_wrapup(2)
add_stim_bar(2, stim_frs, [0, 0, 0])

figure(3)
plot(PMT_traces)
set_xlabels_time(3, 0.03, 10);
fig_wrapup(3)
add_stim_bar(3, stim_frs, [0, 0, 0])

figure(4)
plot(dff_PMT_traces(:, 1:3))
set_xlabels_time(4, 0.03, 10);
fig_wrapup(4)
add_stim_bar(4, stim_frs, [0, 0, 0])

std(dff_MPPC_traces(1:500, 1:3))

std(dff_PMT_traces(1:500, 1:3))


