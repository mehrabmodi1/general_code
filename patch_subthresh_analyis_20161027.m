clear all
close all

single_tr_plot = 0;

direc = 'D:\Data\Janelia\Patch\Data_MM\thacq_files\';
%list_file = 'cell_list_ABs.xls';
%list_file = 'cell_list_ApBp.xls';
%list_file = 'cell_list_G.xls';
%list_file = 'cell_list_unknown.xls';
%list_file = 'cell_list_unknown_unstained.xls';
list_file = 'cell_list_allKCs.xls';

odor_list = {'3-Octanol', ...
             '1-Hexanol', ...
             'Pentyl acetate', ...
             '4-Methylcyclohexanol', ...
             '2-Heptanone', ...
             'Diethyl succinate', ...
             'Ethyl lactate', ...
             '1-Octen-3-ol', ...
             'Geranyl acetate', ...
             'Ethyl acetate', ...
             'Empty' ...
            };



[del del1 cell_list] = xlsread([direc list_file]);

n_cells = size(cell_list, 1);
n_odors = size(odor_list, 2);

v_trace_f_mat = [];
saved_v_areas = [];


for cell_n = 1:n_cells
    cell_path = cell_list{cell_n, 1};
    cell_path = [cell_path, '\'];
    cell_name = cell_list{cell_n, 2};
    
    try

        cell_data = load([cell_path, cell_name, '.mat']);
        spike_data = load([cell_path, cell_name, '_spike.mat']);
        bad_tr_list = load([cell_path, cell_name, '_selection.mat']);
        bad_tr_list = bad_tr_list.rejectedsweeps;
        sf = cell_data.Data_1.parameter.ai_sr;          %AI sampling rate
        s_time = (1./sf);                         %sample time in ms
    catch
        continue
    end
    n_trials = length(fieldnames(cell_data)); 
    
    %initialising stim table, data table
    stim_mat = zeros(n_trials, 3) + nan;        %trial_n x [odor number, odor duration, odor onset time]
    
    
    %loop to load each trial into memory
    for trial_n = 1:n_trials
        
        %skipping trials rejected in THview
        if bad_tr_list(trial_n) == 1
            continue
        else
        end
            
        eval(['curr_data = cell_data.Data_' int2str(trial_n) ';']);
        
        %adding to stim table
        for odor_ni = 1:n_odors
            list_name = odor_list{1, odor_ni};
            curr_name = curr_data.odor;
            if strcmp(list_name, curr_name) == 1
                curr_odor_n = odor_ni;              %current odor number
                stim_mat(trial_n, 1) = curr_odor_n;
                break
            else
            
            end
        end
        
        od_dur = curr_data.parameter.odorD;       %odor duration
        od_onset_t = curr_data.parameter.preO;    %odor onset time
        acq_dur = curr_data.paratable{2, 3};      %trial duration
        
        %skipping short odor presentation trials
        if od_dur == 1
            continue
        else
        end
        
        stim_mat(trial_n, 2:4) = [od_dur, od_onset_t, acq_dur];
        
        %Getting V waveform
        v_trace = curr_data.data.voltage*1000/curr_data.amplifier.Vgain;
        
        %subtracing baseline from V-trace
        base_trace = v_trace(1:(od_onset_t./s_time));
        base_m = mean(base_trace);
        
        %ignoring trials where baseline Vmem was > -50 mV, indicating
        %unhealthy cells
        if base_m > -50
            continue
            
        else
        end
        
        v_trace_orig = v_trace;
        v_trace = v_trace - base_m;
        
        
        %low-pass filtering trace to get rid of spike waveforms
        disp(['trial ' int2str(trial_n) ' of ' int2str(n_trials)])
        v_trace_ff = tsmovavg_m(v_trace,'s', round(.2./s_time), 1);
        v_trace_orig_ff = tsmovavg_m(v_trace_orig, 's', round(.2./s_time), 1);
        
        if single_tr_plot == 1
            
            figure(1)
            plot(v_trace)
            add_stim_shading(1, [od_onset_t./s_time, (od_onset_t./s_time + od_dur./s_time)], 0.4, [.65, .65, .65]);

            figure(2)
            t_vec = (1./sf):(1./sf):acq_dur;
            plot(t_vec, v_trace, 'LineWidth', 2)
            hold on
            %plot(v_trace_f, 'r')
            plot(t_vec, v_trace_ff, 'g', 'LineWidth', 2)
            xlabel('time (s)')
            ylabel('membrane voltage (mV)')
            hold off
        else
        end
        
        %saving filtered versions of all 60s traces
        v_trace_f_mat = [v_trace_f_mat, v_trace_ff];
        
        %saving filtered, background subtracted areas under curve for sus,
        %on time windows for each trial
        on_win = [(od_onset_t./s_time), ((od_onset_t + 2)./s_time)];
        sus_win = [(od_onset_t./s_time), (od_onset_t./s_time + od_dur./s_time)];
        
        a_on = mean(v_trace_orig_ff(on_win(1):on_win(2) ) );          %area under curve for on period
        a_sus = mean(v_trace_orig_ff(sus_win(1):sus_win(2) ) );       %area under curve for sus period
        a_base = nanmean(v_trace_orig_ff(1: (on_win(1) - 1) ) );
        
        %statistical testing
        [h_on, p_on] = ttest2(v_trace_orig_ff(1: (on_win(1) - 1) ), v_trace_orig_ff(on_win(1):on_win(2) ) );
        [h_sus, p_sus] = ttest2(v_trace_orig_ff(1: (on_win(1) - 1) ), v_trace_orig_ff(sus_win(1):sus_win(2) ) );
        
        saved_v_areas = [saved_v_areas; a_base, a_on, a_sus, h_on, h_sus, p_on, p_sus];
        
    end
end

%PLOTTING

%plotting mean filtered, response trace - analogous to Stopfer plot
figure(1)
ftrace_m = nanmean(v_trace_f_mat, 2);
ftrace_se = nanstd(v_trace_f_mat, [], 2)./sqrt(size(v_trace_f_mat, 2));
t_vec = (1./sf):(1./sf):acq_dur;
shadedErrorBar(t_vec,ftrace_m,ftrace_se,{'g'},1)
xlabel('time (s)')
ylabel('membrane voltage (mV)')

figure(4)
plot(t_vec, v_trace_f_mat, 'Color', [0.65, 0.65, 0.65])
hold on
plot(t_vec, v_trace_f_mat(:, 8), 'Color', 'r')
xlabel('time (s)')
ylabel('membrane voltage (mV)')
hold off

%plotting saved v-areas under curve for on vs base periods
figure(2)
plot(saved_v_areas(:, 1), saved_v_areas(:, 2), 'O', 'markerfacecolor', 'b')
xlabel('mean baseline voltage (mV)')
ylabel('mean on-period voltage (mV)')
hold on
sig_cells = find(saved_v_areas(:, 4) == 1);
plot(saved_v_areas(:, 1), saved_v_areas(:, 2), 'O', 'markerfacecolor', 'r')
hold off
ax = axis;
a = min(ax(1, [1, 3]));
b = max(ax(1, [2, 4]));
axis([a, b, a, b]);


%plotting saved v-areas under curve for on vs base periods
figure(3)
plot(saved_v_areas(:, 1), saved_v_areas(:, 3), 'O', 'markerfacecolor', 'b')
xlabel('mean baseline voltage (mV)')
ylabel('mean sustained-period voltage (mV)')
hold on
sig_cells = find(saved_v_areas(:, 5) == 1);
plot(saved_v_areas(:, 1), saved_v_areas(:, 3), 'O', 'markerfacecolor', 'r')
hold off
ax = axis;
a = min(ax(1, [1, 3]));
b = max(ax(1, [2, 4]));
axis([a, b, a, b]);