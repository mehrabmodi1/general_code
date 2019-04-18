clear all
close all

direc = 'D:\Data\CSHL\20160304\analog_acq\';
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');


data_mat = load([direc, 'acqn_2.mat']);


data_mat = data_mat.data_mat;

n_odors = 2;

%cross correlating each trigger trace with all the others to line up all
%the PID traces with the same lags.
trig_trace_mat = squeeze(data_mat(:, :, 2));
pid_trace_mat = squeeze(data_mat(:, :, 1));

n_points = size(data_mat, 1);


%taking every odd numbered trace as odor1 traces and even numbered trace as
%odor2 traces
odor1_trig_traces = trig_trace_mat(:, 1:2:size(trig_trace_mat, 2));
odor1_pid_traces = pid_trace_mat(:, 1:2:size(pid_trace_mat, 2));

odor2_trig_traces = trig_trace_mat(:, 2:2:size(trig_trace_mat, 2));
odor2_pid_traces = pid_trace_mat(:, 2:2:size(pid_trace_mat, 2));

trig_trace_mat = zeros(size(odor1_trig_traces, 1), size(odor1_trig_traces, 2), 2);
pid_trace_mat = zeros(size(odor1_trig_traces, 1), size(odor1_trig_traces, 2), 2);

trig_trace_mat(:, :, 1) = odor1_trig_traces;
trig_trace_mat(:, :, 2) = odor2_trig_traces;

pid_trace_mat(:, :, 1) = odor1_pid_traces;
pid_trace_mat(:, :, 2) = odor2_pid_traces;

n_traces = size(pid_trace_mat, 2);

trig_sorted = zeros(size(trig_trace_mat, 1), n_traces, n_odors) + nan;
pid_sorted = zeros(size(trig_trace_mat, 1), n_traces, n_odors) + nan;

for odor_n = 1:n_odors
    trig_traces = squeeze(trig_trace_mat(:, :, odor_n));
    pid_traces = squeeze(pid_trace_mat(:, :, odor_n));
    lags = xcorr(trig_traces);

    n_traces = size(trig_traces, 2);
    maxlag_mat = zeros(n_traces, n_traces) + nan;             %matrix of lags with maximum correlation
    for trace_ni = 1:n_traces
        for trace_nj = 1:n_traces
            a = lags(:, ((trace_ni-1).*6)+trace_nj );
            [del, maxi] =  max(a);
            maxlag_mat(trace_ni, trace_nj) = maxi;
            
        end
    end
    

    %applying lags to trigger traces and odor traces and re-assembling matrices
    %for both
    
    for trace_n = 1:n_traces
        cur_trig_trace = trig_traces(:, trace_n);
        cur_pid_trace = pid_traces(:, trace_n);
        tr_end = length(cur_trig_trace);
        cur_lag = maxlag_mat(1, trace_n) - tr_end;
        lag_sign = sign(cur_lag);
        cur_lag = abs(cur_lag);
        if lag_sign == 1
            trace_insert_trig = cur_trig_trace(1:(tr_end - cur_lag + 1));
            trace_insert_pid = cur_pid_trace(1:(tr_end - cur_lag + 1));
            trig_sorted( (cur_lag + 1):tr_end, trace_n, odor_n) = trace_insert_trig;
            pid_sorted( (cur_lag + 1):tr_end, trace_n, odor_n) = trace_insert_pid;
            
        elseif lag_sign == 0
            trace_insert_trig = cur_trig_trace( (cur_lag + 1):tr_end);
            trace_insert_pid = cur_pid_trace( (cur_lag + 1):tr_end);
            trig_sorted( (cur_lag + 1):tr_end, trace_n, odor_n) = trace_insert_trig;
            pid_sorted( (cur_lag + 1):tr_end, trace_n, odor_n) = trace_insert_pid;
        elseif lag_sign == -1
            trace_insert_trig = cur_trig_trace(cur_lag:tr_end);
            trace_insert_pid = cur_pid_trace(cur_lag:tr_end);
            trig_sorted(1:(tr_end - cur_lag + 1), trace_n, odor_n) = trace_insert_trig;
            pid_sorted(1:(tr_end - cur_lag + 1), trace_n, odor_n) = trace_insert_pid;
            
        end
        
        

    end
    
    
end

figure(1)
plot(trig_sorted(:, :, 1))
hold on
plot(pid_sorted(:, :, 1), 'LineWidth', 2, 'Color', color_vec(6, :))
xlabel('time')
ylabel('PID signal (V)')
axis([0, 60000, 0, 1])
set_xlabels_time(1, 10e-4, 0.5)

figure(2)
plot(trig_sorted(:, :, 2))
hold on
plot(pid_sorted(:, :, 2), 'LineWidth', 2, 'Color', color_vec(3, :))
xlabel('time')
ylabel('PID signal (V)')
set_xlabels_time(2, 10e-4, 0.5)
%axis([0, 60000, 0, 0.8])
