function [resp_sizes, sig_trace_mat, sig_cell_mat, sig_cell_mat_key, resp_areaundercurves] = cal_sig_responses_res_modular(dff_data_mat, stim_mat_struc, stim_mat_simple, frame_time, od_col_ns, dur_col_ns)
%This function takes as inputs, the 3-D dff_data_mat that contains dF/F traces stored
%by frame_n, cell_n, trial_n and odor_n; and stim_mat which contains
%stimulus delivery information for each trial. The outputs are three 2D
%matrices. Resp_areas of size n_cell, n_trials contains the area under the pk response 
%of each cell on each trial sig_trace_mat of size n_cells, n_trials indicates which
%individual traces were significant responses, and sig_cell_mat of size 
%n_cells, n_odors which indicates which cells responded for more than half
%the presentations of any single odor.

%replacing nans in stim_mat_simple with 0s
del = isnan(stim_mat_simple);
stim_mat_simple(del) = 0;

del = find(stim_mat_simple(:, dur_col_ns(1)) < 1);
stim_mat_simple(del, dur_col_ns(1)) = 0;

n_frames = size(dff_data_mat, 1);
n_cells = size(dff_data_mat, 2);
n_trials = size(dff_data_mat, 3);
od_list_olf1 = unique(stim_mat_simple(:, od_col_ns(1) ));
od_list_olf2 = unique(stim_mat_simple(:, od_col_ns(2) ));
dur_list_olf1 = unique(stim_mat_simple(:, dur_col_ns(1) ));
dur_list_olf2 = unique(stim_mat_simple(:, dur_col_ns(2) ));
dur_list_olf1 = unique([0; dur_list_olf1]);
od_list_olf2 = unique([0; od_list_olf2]);
dur_list_olf2 = unique([0; dur_list_olf2]);

pre_stim_frame = 50;

resp_sizes = zeros(n_cells, n_trials) + nan;
resp_areaundercurves = zeros(n_cells, n_trials) + nan;
sig_trace_mat = zeros(n_cells, n_trials);

gp_list = [];
%builing list of odor, odor-dur pairs

for olf1_dur_n = 1:length(dur_list_olf1)
    olf1_dur_ni = dur_list_olf1(olf1_dur_n);
    
    for olf2_dur_n = 1:length(dur_list_olf2)
            olf2_dur_ni = dur_list_olf2(olf2_dur_n);
        
        for olf1_od_n = 1:length(od_list_olf1)
            olf1_od_ni = od_list_olf1(olf1_od_n);
        
            for olf2_od_n = 1:length(od_list_olf2)
                olf2_od_ni = od_list_olf2(olf2_od_n);
             
                gp_list = [gp_list; olf1_od_ni, olf1_dur_ni, olf2_od_ni, olf2_dur_ni];
            end
        end
        
    end
end 

sig_cell_mat = [];
sig_cell_mat_key = [];
for rep_gp = 1:size(gp_list, 1)               %a rep_gp is a group of repeats of the same stimulus types
    curr_gp = gp_list(rep_gp, :);       %stimuli delivered in current group of repeats
    
    olf1_od_ni = curr_gp(1);
    olf1_dur_ni = curr_gp(2);
    olf2_od_ni = curr_gp(3);
    olf2_dur_ni = curr_gp(4);
    
    rep_tr_list = find(stim_mat_simple(:, od_col_ns(1) ) == olf1_od_ni & stim_mat_simple(:, dur_col_ns(1)) == olf1_dur_ni &...
        stim_mat_simple(:, od_col_ns(2)) == olf2_od_ni & stim_mat_simple(:, dur_col_ns(2)) == olf2_dur_ni);    %list of repeat tr numbers
    
    if isempty(rep_tr_list) == 1
        continue
    else
    end
    sig_cell_mat_key = [sig_cell_mat_key; curr_gp];
    
    %identifying odor stimulus periods for olf1 and olf2, and including
    %beginning of first pulse to end of last pulse as analysis period.
    stim_frs = compute_stim_frs_modular(stim_mat_struc, rep_tr_list(1), frame_time);
    stim_frs = [stim_frs{1}, stim_frs{2}];
    stim_frs = [min(stim_frs, [], 'omitnan'), max(stim_frs, [], 'omitnan')];
    sig_cell_vec = zeros(n_cells, 1);    
    for cell_n = 1:n_cells
        curr_traces = squeeze(dff_data_mat(:, cell_n, rep_tr_list));             
        
        %moving window filtering each trace
        curr_traces_f = zeros(size(curr_traces, 1), size(curr_traces, 2)) + nan;
        for rep_n = 1:size(curr_traces, 2)
            curr_trace = curr_traces(:, rep_n);
            curr_traces_f(:, rep_n) = movmean(curr_trace, round(0.5./frame_time));        %filtering with a 0.5 s wide box-car
        end
               
        base_traces = curr_traces_f(pre_stim_frame:(stim_frs(1) - 2), :);                 %pre-odor-stim baseline frames
        base_m = median(base_traces, 1, 'omitnan');
        base_s = std(base_traces, 1, 'omitnan');
        threshes = abs(base_m) + 2.33.*base_s;
        thresh_mean = mean(mean(base_traces, 2, 'omitnan'), 'omitnan') + 2.33.*std(mean(base_traces, 2, 'omitnan'), 'omitnan');       %significance criterion for resp in mean trace
        
        try
            resp_traces = curr_traces_f(stim_frs(1):round(stim_frs(2) + 5./frame_time), :); %extending analysis window to 5s after odor off to capture off responses
        catch
            keyboard
        end
        
        ave_trace = mean(resp_traces, 2, 'omitnan');
        
        %identifying peak in rep-averaged trace to assign time of pk-response
        %this enforces a fixed pk time for each response in the smoothed trace
        [ave_pk, ave_pk_fr] = max(ave_trace, [], 'omitnan');
        %setting criterion that ave trace should also have a sig response
        if ave_pk < thresh_mean
            continue
        else
            
        end
        
        pk_resps = resp_traces(ave_pk_fr, :);                  %single repeat response amplitudes at pk time of avg trace
        %resp_areas(cell_n, rep_tr_list) = pk_resps;
        try
            resp_traces_win = resp_traces((max([(ave_pk_fr-round(1./frame_time)), 1])):(min([(ave_pk_fr+round(1./frame_time)), size(resp_traces, 1)])), :);
            resp_sizes(cell_n, rep_tr_list) = max(resp_traces_win, [], 1);
            resp_areaundercurves(cell_n, rep_tr_list) = mean(resp_traces, 1, 'omitnan');
        catch
            keyboard
        end
        pk_resps_indiv = max(resp_traces);
        
        
        %checking if each response is significantly higher than baseline
        for rep_n = 1:size(pk_resps, 2) 
            curr_pk = pk_resps(rep_n);
            if curr_pk > threshes(rep_n)
                
                sig_trace_mat(cell_n, rep_tr_list(rep_n)) = 1;
                
%                 %visualising current trace to check if its a sig resp
%                 figure(1)
%                 curr_trace = curr_traces_f(:, rep_n);
%                 plot(curr_trace)
%                 hold on
%                 thresh_vec = repmat(threshes(rep_n), size(curr_traces, 1), 1);
%                 plot(thresh_vec, 'r')
%                 plot(base_traces(:, rep_n), 'g')
%                 plot((ave_pk_fr + length(base_traces)), curr_pk, 'Or')
%                 ave_trace = nanmean(curr_traces_f, 2);
%                 plot(ave_trace, 'k', 'LineWidth', 2)
%                 hold off
%                 keyboard
%                 
            else
                
            end
            
            
            %in parallel checking if single trace peaks (not values at avg
            %trace pk) are significant as per the old criterion for 1s
            %stims
            if pk_resps_indiv(rep_n) > threshes(rep_n)
                sig_trace_mat_old(cell_n, rep_tr_list(rep_n)) = 1;
            else
                
            end
                
        end
        
        %deciding if the current cell is a significant responder in the current rep-gp
        %if sum(sig_trace_mat(cell_n, rep_tr_list))./length(rep_tr_list) >= 0.5      %at least half the repeats should be sigonificant resps
            
        if sum(sig_trace_mat(cell_n, rep_tr_list)) > 2                          %at least 2 of the repeats shold be significant resps
            sig_cell_vec(cell_n, 1) = 1;
            h = 1;             
        else
            h = 0;
        end
        
        
        
%         if h == 1
%             n_points = sum(abs(isnan(curr_traces_f(:, 1)) - 1));                     %not plotting padding nans
%             if n_points == 0
%                 n_points = sum(abs(isnan(curr_traces_f(:, 2)) - 1));                     %not plotting padding nans
%             else
%             end
%             thresh_vec = repmat(threshes, n_points, length(threshes));
%             
%             figure(2)
%             plot(curr_traces_f(1:n_points, :), 'Color', [0.65, 0.65, 0.65])
%             hold on
%             ave_full_trace = nanmean(curr_traces_f(1:n_points, :), 2);
%             plot(ave_full_trace, 'Color', 'k', 'LineWidth', 2)
%             plot(thresh_vec, 'r')
%             add_stim_shading(2, [stim_frame, stim_end_fr], 0.4, [.65, .65, .65]);
%             hold off
%             set_xlabels_time(2, frame_time./1000, 1.5)
%             figure(2)
%             del = input('press enter for next sig cell');
%             hold off
%         else
%         end
    end
    sig_cell_mat = [sig_cell_mat, sig_cell_vec];    
    
end
