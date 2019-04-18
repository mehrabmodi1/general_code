clear all
close all

direc = 'D:\Data\CSHL\20160622\analog_acq\acq3\';
dir_contents = dir(direc);

aqn_rate = 1000;       %DAQ acqn rate in Hz

trace_file = dir_contents(3).name;
param_file = dir_contents(4).name;

trace_mat = load([direc, trace_file]);
trace_mat = trace_mat.data_mat;

param_mat = load([direc, param_file]);
param_mat = param_mat.params;


odor_trials = param_mat.odours;     %odor delivered on each trial
odor_list = unique(odor_trials);    %list of odors delivered

n_trials = length(odor_trials);
n_repeats = param_mat.reps;

n_pulses_trials = param_mat.n_od_pulses;
n_pulses_list = unique(n_pulses_trials);
ipi_trials = param_mat.inter_pulse_interval;
ipi_list = unique(ipi_trials);
duration_trials = param_mat.duration;
duration_list = unique(duration_trials);

%checking if product of all parameter ranges equals number of trials 
if (length(n_pulses_list).*length(odor_list).*length(ipi_list).*length(duration_list)) ~= (n_trials./n_repeats)
    disp('Warning, parameter ranges dont square off against number of trials.')
else
end

saved_info = [];
for odor_n = 1:length(odor_list)
    odor_ni = odor_list(odor_n);
    od_tr_list = find(odor_trials == odor_ni); 
    for n_pulses_n = 1:length(n_pulses_list)
        n_pulses_ni = n_pulses_list(n_pulses_n);
        n_pulses_tr_list = find(n_pulses_trials == n_pulses_ni);
        for ipi_n = 1:length(ipi_list)
            ipi_ni = ipi_list(ipi_n);
            ipi_tr_list = find(ipi_trials == ipi_ni);
            for duration_n = 1:length(duration_list)
                duration_ni = duration_list(duration_n);
                duration_tr_list = find(duration_trials == duration_ni);
                
                curr_tr_list = intersect(od_tr_list, n_pulses_tr_list);
                curr_tr_list = intersect(curr_tr_list, ipi_tr_list);
                curr_tr_list = intersect(curr_tr_list, duration_tr_list);
                
                stim_onset = param_mat.stimLatency(curr_tr_list(1)).*aqn_rate;
                stim_end = stim_onset + ( (duration_ni + ipi_ni).*n_pulses_ni - ipi_ni ).*aqn_rate;
               
                curr_traces = trace_mat(:, curr_tr_list);
                
                %baseline subtraction
                baseline_vec = mean(curr_traces(1:(stim_onset-1), :), 1);
                curr_traces = curr_traces - repmat(baseline_vec, size(curr_traces, 1), 1);
                
                area_vec = sum(curr_traces(stim_onset:stim_end, :), 1);
                area_mean = mean(area_vec);
                area_sd = std(area_vec);
                sd_percentage = (area_sd./area_mean).*100;
                
                figure(1)
                plot([(1./aqn_rate):(1./aqn_rate):(size(curr_traces, 1)./aqn_rate)], curr_traces)
                xlabel('time (s)')
                ylabel('PID signal (AU)')
                title(['odor ' int2str(odor_ni) '; duration ' num2str(duration_ni), '; n pulses ' num2str(n_pulses_ni) '; ipi ' num2str(ipi_ni) '; area SD ' num2str(sd_percentage) '%' ])
                del = input('press enter')
                
                saved_info = [saved_info; area_mean, area_sd, sd_percentage];
                
            end
        
        end
        
    end
        
end

