clear all
close all

%loading in PID traces for alchohol concentration series
alch_path = 'C:\Data\Data\Raw_data\20180605\olf_testing\';
[PID_traces] = get_PID_traces(alch_path, [1:50], 1./2000);
PID_traces = movmean(PID_traces', 0.05.*2000)';     %moving window averaging traces with a 50 ms window

[stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains(alch_path, []);

conc_list = unique(stim_mat_simple(:, 8));
fin_conc_list = conc_list.*0.0204;
mean_hts = zeros(5, 1) + nan;
for conc_n = 1:5
    conc_ni = conc_list(conc_n);
    curr_conc = conc_ni.*(0.0204);
    curr_trs = find(stim_mat_simple(:, 8) == conc_ni);
    ave_trace = mean(PID_traces(curr_trs, :), 1, 'omitnan');
    mean_hts(conc_n) = mean(ave_trace(5.618e4:6.606e4));
    
    figure(1)
    plot(ave_trace, 'lineWidth', 3, 'Color', ([0.9, 0.9, 0.9].*(1 - conc_ni)))
    hold on
    
    legend('1.6%', '1.2%', '0.8%', '0.4%', '0.2%', 'northwest')
    ylabel('Odor concentration (A.U.)')
    set_xlabels_time(1, 1./2000, 10)
end

%plotting linearity fit curve
fit = polyfit(fin_conc_list, mean_hts, 1);
fit_pts = polyval(fit, [0; fin_conc_list]);

figure(2)
plot((fin_conc_list.*100), mean_hts, 'r.', 'markerSize', 14)
hold on
plot([0; (fin_conc_list.*100)], fit_pts, 'b--')
xlabel('deliv. odor conc. (% of saturation)')
ylabel('meas. odor conc. (AU)')
axis([0, 2, 0, 0.9])
fig_wrapup(3)

%loading in PID traces for PA BA and EL at 0.2% concentration
od_path = 'C:\Data\Data\Raw_data\20180605\olf_testing_PaBaEl\';
[PID_traces] = get_PID_traces(od_path, [1:30], 1./2000);
[stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains(od_path, []);
PID_traces = movmean(PID_traces', 0.05.*2000)';         %moving window averaging traces with a 50 ms window
odor_list = unique(stim_mat_simple(:, 2));

for odor_n = 1:3
   odor_ni = odor_list(odor_n);
   curr_trs = find(stim_mat_simple(:, 2) == odor_ni);
   
   ave_trace = mean(PID_traces(curr_trs, :), 1, 'omitnan');
   sd_trace = std(PID_traces(curr_trs, :), 1, 'omitnan');
   
   figure((2 + odor_n))
   %plot(ave_trace, 'lineWidth', 2, 'Color', color_vec(odor_n, :))
   %errorbar(1:1:length(ave_trace), ave_trace, sd_trace, 'Color', color_vec(odor_n, :))
   shadedErrorBar(1:1:length(ave_trace),ave_trace,sd_trace,{'Color', color_vec(odor_n, :)}, 0)
   
   ylabel('odor concentration (A.U.)')
   set_xlabels_time((2 + odor_n), 1./2000, 10)
   fig_wrapup((2+odor_n))
end

fig_wrapup(1)
fig_wrapup(2)