clear all
close all

direc = 'D:\Data\Janelia\resonant\20171213\olf_calib3';

color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');

param_mat = load([direc, '\params.mat']);
param_mat = param_mat.params_mat;

[stim_mat_simple, column_heads, rand_trains] = load_params_res2(param_mat);
n_trials = size(stim_mat_simple, 1);
odor_list = sort(unique(stim_mat_simple(:, 1)));
odor_list(3) = [];

n_odors = length(unique(odor_list));
n_rand_trains = max(stim_mat_simple(:, 10));
del = find(stim_mat_simple(:, 1) == 1);
n_reps = length(del);
stim_point = stim_mat_simple(1, 6).*2000;        %PID datapoint number when stim begins. stim_latency multiplied by the acq rate ie. 2KHz
stim_end_point = (stim_mat_simple(1, 2).*2000) + stim_point;
n_trains = max(stim_mat_simple(:, 10));

mean_trace_mat = [];
sd_trace_mat = [];
for odor_n = 1:n_odors
    odor_ni = odor_list(odor_n);
    
    for rand_train_n = 1:n_rand_trains
        curr_trs = find(stim_mat_simple(:, 1) == odor_ni & stim_mat_simple(:, 10) == rand_train_n);
        trace_mat = [];
        %loading in traces for current odor
        for tr_n = 1:length(curr_trs)
            tr_ni = curr_trs(tr_n);
            trace = load([direc, '\PID_trace_tr-', int2str(tr_ni)]);
            trace = trace.PID_data;

            if tr_n == 1
                trace_mat = [trace_mat, trace(:, 1)];
                time_vec = trace(:, 2);
            elseif tr_n > 1
                trace_mat = pad_n_concatenate(trace_mat, trace(:, 1), 2, nan);
            end
        end
    
        %calculating dF/F
        mean_baselines = mean(trace_mat(1:stim_point, :), 1, 'omitnan');
        mean_baselines = repmat(mean_baselines, size(trace_mat, 1), 1);
        trace_mat = (trace_mat - mean_baselines);
        
        trace_mat_all(:, :, rand_train_n, odor_n) = trace_mat;
        mean_trace_mat = [mean_trace_mat, mean(trace_mat, 2, 'omitnan')];
        sd_trace_mat = [sd_trace_mat, std(trace_mat, 0, 2, 'omitnan')];
        
    end
    

end

%Analysing traces
%calculating conc dips over stim train
init_pks = max(mean_trace_mat(stim_point:(stim_point + 10000), :));         %conc early in stim train
late_pks = max(mean_trace_mat((stim_end_point - 10000):stim_end_point, :)); %conc late in stim train
percentage_drops = (init_pks - late_pks)./init_pks;                         %drop off in conc from early in stim train to end of stim train

%identifying shortest pulses
for train_n = 1:n_trains
    tr_n = find(stim_mat_simple(:, 10) == train_n);
    tr_n = tr_n(1);
    stim_train = param_mat(tr_n).rand_train;
    s_latency = stim_mat_simple(1, 6);
    
    %sampling points around shortest odor pulse
    [del, short_pulse] = min(stim_train(:, 2));
    if short_pulse > 1
        pulse_on_t = stim_train(short_pulse, 1) + s_latency + sum(sum(stim_train(1:(short_pulse - 1), :)));
    else
        pulse_on_t = stim_train(short_pulse, 1) + s_latency;
    end
    pulse_on_pt = pulse_on_t.*2000;
    
    trace_n_vec = train_n:n_trains:(n_odors.*n_trains);                   %vector of mean trace numbers for current odor train across all odors
    del = mean_trace_mat( (pulse_on_pt - 4000):(pulse_on_pt + 16000), trace_n_vec);
    short_mean_traces = zeros(length(del), size(del, 2)) + nan;
    short_mean_traces(1:7000, :) = del(1:7000, :); 
    del = sd_trace_mat( (pulse_on_pt - 4000):(pulse_on_pt + 16000), trace_n_vec);
    short_sd_traces = zeros(length(del), size(del, 2)) + nan;
    short_sd_traces(1:7000, :) = del(1:7000, :);
    
    %sampling points around longest odor pulse
    [del, long_pulse] = max(stim_train(:, 2));
    if long_pulse > 1
        pulse_on_t = stim_train(long_pulse, 1) + s_latency + sum(sum(stim_train(1:(long_pulse - 1), :)));
    else
        pulse_on_t = stim_train(long_pulse, 1) + s_latency;
    end
    pulse_on_pt = pulse_on_t.*2000;
    
    long_mean_traces = mean_trace_mat( (pulse_on_pt - 4000):(pulse_on_pt + 16000), trace_n_vec);       
    long_sd_traces = sd_trace_mat( (pulse_on_pt - 4000):(pulse_on_pt + 16000), trace_n_vec);
    
    %plotting
    %plotting entire pulse trains
    figure(1)
    t_vec = -1.*s_latency:0.0005:(0.0005.*size(mean_trace_mat, 1) - s_latency - 0.0005);
    for plot_n = 1:size(short_mean_traces, 2)
        shadedErrorBar(t_vec, mean_trace_mat(:, trace_n_vec(plot_n)), sd_trace_mat(:, trace_n_vec(plot_n)), {'Color', color_vec(plot_n, :)}, 0)
%         errorbar(t_vec, mean_trace_mat(:, trace_n_vec(plot_n)), sd_trace_mat(:, trace_n_vec(plot_n)), 'Color', color_vec(plot_n, :))
        
        hold on
    end
    drawnow
    hold off
    
    %zoomed in on shortest, longest pulses
    figure(2)
    t_vec = -2:0.0005:(0.0005.*size(short_mean_traces, 1) - 2);
    t_vec = t_vec(2:end);
    
    for plot_n = 1:size(short_mean_traces, 2)
        shadedErrorBar(t_vec', short_mean_traces(:, plot_n), short_sd_traces(:, plot_n), {'Color', color_vec(plot_n, :)}, 0)
        hold on
        shadedErrorBar(t_vec', long_mean_traces(:, plot_n), long_sd_traces(:, plot_n), {'Color', color_vec(plot_n, :)}, 0)
         
%         errorbar(t_vec', short_mean_traces(:, plot_n), short_sd_traces(:, plot_n), 'Color', color_vec(plot_n, :))
%         hold on
%         errorbar(t_vec', long_mean_traces(:, plot_n), long_sd_traces(:, plot_n), 'Color', color_vec(plot_n, :))
    end
    drawnow
    hold off
    keyboard
end

