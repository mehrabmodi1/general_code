clear all
close all

direc = 'D:\Data\Janelia\Patch\Data_MM\thacq_files\';
%list_file = 'cell_list_ABs.xls';
%list_file = 'cell_list_ApBp.xls';
%list_file = 'cell_list_G.xls';
list_file = 'cell_list_unknown.xls';
%list_file = 'cell_list_unknown_unstained.xls';

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

plot_spike_trace = 0;
plot_single_cell_stuff = 0;

saved_rates = [];
saved_PSTH_curves = [];
saved_sig_cells = [];
saved_sig_cells2 = [];

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
    catch
        continue
    end
    n_trials = length(fieldnames(cell_data)); 
    
    %initialising stim table, data table
    stim_mat = zeros(n_trials, 3) + nan;        %trial_n x [odor number, odor duration, odor onset time]
    sp_data_mat = zeros(800000, n_trials) + nan;
    sp_wav_mat = [];        
    %loop to load each trial into memory
    for trial_n = 1:n_trials
        
        %skipping trials rejected in THview
        if bad_tr_list(trial_n) == 1
            continue
        else
        end
            
        eval(['curr_data = cell_data.Data_' int2str(trial_n) ';']);
        
        try
            eval(['sp_times = spike_data.Spike_' int2str(trial_n) ';']);
            sp_datapoints = sp_times(:, 1).*sf;
        catch
            sp_times = [];
            sp_datapoints = [];
        end
        
        
        %sampling spike waveforms
        v_trace = curr_data.data.voltage*1000/curr_data.amplifier.Vgain;
        if plot_spike_trace == 1
            warning('off', 'MATLAB:colon:nonIntegerIndex');
                    for sp_n = 1:size(sp_datapoints, 1);
                sp_point = sp_datapoints(sp_n);
                sp_wv = v_trace(max((sp_point - (.05.*sf) ), 1):min((sp_point + (.05.*sf) ), length(v_trace) ));
                sp_wv = sp_wv - mean(sp_wv);
                sp_wav_mat = pad_n_concatenate(sp_wav_mat, (sp_wv)', 1, nan);
                    end
            warning('on', 'MATLAB:colon:nonIntegerIndex');
        else
        end
        
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
        
        stim_mat(trial_n, 2) = curr_data.parameter.odorD;       %odor duration
        stim_mat(trial_n, 3) = curr_data.parameter.preO;        %odor onset time
        stim_mat(trial_n, 4) = curr_data.paratable{2, 3};       %trial duration
        
        v_trace = curr_data.data.voltage;
        sp_vec = zeros((curr_data.parameter.dur.*sf), 1);
        try
            sp_vec(round(sp_datapoints), 1) = 1;
        catch
            keyboard
        end
        sp_data_mat(1:length(sp_vec), trial_n) = sp_vec;
    end
    n_total_spikes = nansum(nansum(sp_data_mat)); 
    
        
    %plotting all spike waveforms for current cell
    if plot_spike_trace == 1
        base_color = [.4, .5, .9];
        if n_total_spikes > 1
            color_vec = (.05:((1-.05)./400):1)';
            color_vec = repmat(color_vec, 1, 3);
            color_vec(1, :) = [];
            color_vec = color_vec.*repmat(base_color, 400, 1);
        else
            color_vec = base_color;
        end
        fig5 = figure(5);
        t_vec = -50:1./(sf./1000):50;
        for spike_ni = 1:400
            r_sp_list = round(1 + rand(1, 400).*(n_total_spikes - 1) );
            curr_spike = r_sp_list(spike_ni);
            curr_wav = sp_wav_mat(curr_spike, :);
            curr_color = color_vec(spike_ni, :);
            plot(t_vec, curr_wav, 'LineWidth', 1, 'Color', curr_color)
            hold on
        end
        plot(t_vec, nanmean(sp_wav_mat), 'LineWidth', 3, 'Color', [0, 0, 0]);
        hold off
        xlabel('time (ms)');
        ylabel('mean-subtracted voltage (mV)');
        set(fig5, 'Position', [100, 100, 200, 500]);
    else
    end
    
    
    %% Creating raster plots and PSTH curves
    odor_dur_list = unique(stim_mat(:, 2));
    long_dur_list = odor_dur_list(odor_dur_list > 1);               %list of odor durs longer than 1s.
    long_trs = stim_mat(:, 2) > 1;
    long_odor_list = unique(stim_mat(long_trs, 1));            %list of odors delivered for long durations.
    for long_odor_n = 1:length(long_odor_list)
        long_odor = long_odor_list(long_odor_n);
        
        for odor_dur_n = 1:length(odor_dur_list)
            odor_dur = odor_dur_list(odor_dur_n);
            curr_dur_trs = find(stim_mat(:, 2) == odor_dur);      %list of trials of current odor duration.
            curr_trs = find(stim_mat(:, 1) == long_odor);        %list of trials for current odor (which also has long trials).
            curr_trs = intersect(curr_trs, curr_dur_trs);
            
            %skipping if there are no trials - can happen for the odor empty where
            %there are no short trials
            if isempty(curr_trs) == 1
                disp(['skipping odor dur ' int2str(odor_dur)])
                continue
            else
            end
            
            curr_spike_mat = sp_data_mat(:, curr_trs);
            
            %getting rid of padding nans
            trial_dur = stim_mat(curr_trs(1), 4);
            trial_pts = trial_dur.*sf;
            curr_spike_mat((trial_pts + 1):size(curr_spike_mat, 1), :) = [];
                        
            %skipping neuron-odor-duration combinations with 0 spikes
%             if sum(sum(curr_spike_mat)) == 0
%                 disp(['skipping odor dur, no spikes'])
%                 continue
%             else
%             end
            %lag = size(curr_spike_mat, 1)./2000;
            %curr_spike_mat_f = tsmovavg_m(curr_spike_mat,'s',lag,1);      %filtering spikes to spread them out so that each spike is not sub-pixel resolution for display

            if plot_single_cell_stuff == 1
                %creating raster plot
                fig1 = figure(odor_dur_n);
                lag = size(curr_spike_mat, 1)./2000;                          %choosing a lag to make sure there are 2000 lines in the plot.
                curr_spike_mat_f = tsmovavg_m(curr_spike_mat,'s',lag,1);      %filtering spikes to spread them out so that each spike is not sub-pixel resolution for display
                pix_vals = unique(curr_spike_mat_f);
                curr_spike_mat_th = im2bw(curr_spike_mat_f, pix_vals(2)./2);  %re-binarising with a threshold half the lowest value greater than 0
                colormap([1, 1, 1; 0, 0, 0]);
                imagesc(curr_spike_mat_th')
                ylabel('repeats')
                set_xlabels_time(odor_dur_n, 1./sf, .5);


                stim_point = stim_mat(curr_trs(1), 3).*sf;
                stim_end_point = stim_mat(curr_trs(1), 2).*sf + stim_point;
                stim_frs = [stim_point, stim_end_point];
                add_stim_shading(odor_dur_n, stim_frs, .25, [0.6, 0.8, 0.6]);
                set(fig1, 'Position', [100, 100, 500, 125]);
            else
            end

            %Plotting PSTH curve
            tr_dur = size(curr_spike_mat, 1)./sf;        %trial duration in s
            bin_width = 0.5;                               %in s
            n_bins = tr_dur./bin_width;
            bin_width_pts = bin_width.*sf;

            PSTH_curves = zeros(n_bins, size(curr_spike_mat, 2)) + nan;
            for bin_n = 1:n_bins
                bin_vec = [zeros(1, bin_width_pts.*(bin_n-1)), ones(1, bin_width_pts), zeros(1, tr_dur.*sf - (bin_width_pts.*bin_n))];
                bin_count = bin_vec*curr_spike_mat;
                PSTH_curves(bin_n, :) = bin_count;          %PSTH curves are histograms of spike times, convertible to spikes/s in each 500 ms time bin
            end

            PSTH_curves = PSTH_curves.*(1./bin_width);           %converting spike counts from counts/bin to Hz

            if odor_dur == 60
                
                del = mean(PSTH_curves, 2);
                if max(del(40:120)) > 8
                    keyboard
                else
                end
                
                saved_PSTH_curves = [saved_PSTH_curves; mean(PSTH_curves, 2)'];
                
                %testing for sigdiff between sus, pre spike rates
                pre_pts = (PSTH_curves(1:(3./bin_width), :));
                on_pts = (PSTH_curves(((3./bin_width) + 1):(8./bin_width), :));
                sus_pts = (PSTH_curves((8./bin_width + 1):(63./bin_width), :));     %each time bin is an independent sample of the current spike rate
                [h, p] = ttest2(pre_pts, sus_pts);
                %cell is a sus responder if it responds on more than half
                %the trials
                if mean(h) > 0
                    save_vec(1, 1) = 1;     
                else
                    save_vec(1, 1) = 0;
                end
                
                [h, p] = ttest2(pre_pts, on_pts);
                %cell is an on responder if it responds on more than half
                %the trials
                if mean(h) >= 0.5
                    save_vec(1, 2) = 1;     
                else
                    save_vec(1, 2) = 0;
                end
                
                
                saved_sig_cells = [saved_sig_cells; save_vec];
                
            else
            end
                
            if plot_single_cell_stuff == 1
                t_vec = bin_width:bin_width:tr_dur;
                fig3 = figure(odor_dur_n + 2);
                plot(t_vec, PSTH_curves, 'LineWidth', 2, 'Color', [0.9, 0.9, 0.9])
                hold on
                plot(t_vec, nanmean(PSTH_curves, 2), 'LineWidth', 3, 'Color', [0, 0, 0])
                hold off
                xlabel('time (s)')
                ylabel('spike rate (Hz)')
                title(['cell ' int2str(cell_n) ', ' odor_list{long_odor} ', duration ' int2str(odor_dur) ' s' ])

                stim_point = stim_mat(curr_trs(1), 3);
                stim_end_point = stim_mat(curr_trs(1), 2) + stim_point;     %x scale is now s
                stim_pts = [stim_point, stim_end_point];
                add_stim_shading((odor_dur_n + 2), stim_pts, .25, [0.6, 0.8, 0.6]);
                set(fig3, 'Position', [100, 100, 500, 200]);
            else
            end
            
            %keeping track of spike numbers in different time bins and
            %plotting these
            if odor_dur == 60
                pre_bin = [1./sf, (3.*sf - 1)];           %in sample points
                on_bin = [3.*sf, (8.*sf - 1)];            %on period sample points (5s window?)
                sus_bin = [8.*sf, (63.*sf - 1)];          %sus period sample points
                off_bin = [63.*sf, (65.*sf - 1)];         %off period sample points
                
                pre_rate = mean(mean(curr_spike_mat(pre_bin(1):pre_bin(2), :)) );          %mean spikes per second during pre time window
                pre_sd = mean(var(curr_spike_mat(pre_bin(1):pre_bin(2), :)) ).^0.5;        %averaging variances and sqrt-ing to get trial averaged SD 
                on_rate = mean(mean(curr_spike_mat(on_bin(1):on_bin(2), :)) );             %mean spikes per second during on time window
                sus_rate = mean(mean(curr_spike_mat(sus_bin(1):sus_bin(2), :)) );          %mean spikes per second during sus time window
                off_rate = mean(mean(curr_spike_mat(off_bin(1):off_bin(2), :)) );          %mean spikes per second during off time window
                                
                rate_vec = ([pre_rate, on_rate, sus_rate, off_rate] - pre_rate)./pre_sd;   %calculating z-scored spike rates
                saved_rates = [saved_rates; rate_vec];
                
                    
                    
            else
            end
            
        end
            
        %keyboard
    end


   %keyboard 
end

%plotting sus rates versus pre rates
% fig1 = figure(1);
% plot(saved_rates(:, 1), saved_rates(:, 3), 'O', 'MarkerFaceColor', [.4, .6, .9], 'MarkerSize', 8, 'Color', [.4, .6, .9])
% hold on
% max_rate = max([saved_rates(:, 1); saved_rates(:, 3)]);
% plot([0, max_rate], [0, max_rate], '-.', 'Color', [.75, .75, .75], 'LineWidth', 3)
% sig_pts = find(saved_sig_cells(:, 1) == 1);
% %plot(saved_rates(sig_pts, 1), saved_rates(sig_pts, 3), 'O', 'MarkerFaceColor', [.9, .2, .2], 'MarkerSize', 8, 'Color', [.4, .6, .9])
% xlabel('mean spike rate baseline period (Hz)')
% ylabel('mean spike rate sustained period (Hz)')
% hold off
%plotting dist of z-scored sus period rates

fig1 = scattered_dot_plot(saved_rates(:, [2, 3]), 1, 2, 2, 8, [.45, .45, .65], 0, [0, 0, 0], {'on win', 'sus win'});
hold on

ylabel('z-scored spike rate')


%plotting mean PSTH curve (relative of the Stopfer plot)
fig3 = figure(3);
PSTH_mean = nanmean(saved_PSTH_curves);
PSTH_se = std(saved_PSTH_curves)./sqrt(size(PSTH_curves, 2));
shadedErrorBar([0.5:0.5:80], PSTH_mean, PSTH_se, {'Color', '[.4, .6, .9]'});
add_stim_shading(3, [3, 63], .25, [0.6, 0.8, 0.6]);
xlabel('time (s)')
ylabel('average spike rate (Hz)')       %averaged across neuron-odor pairs for 60s duration trials