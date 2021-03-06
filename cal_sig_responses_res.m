function [resp_sizes, sig_trace_mat, sig_trace_mat_old, sig_cell_mat, resp_areaundercurves] = cal_sig_responses_res(dff_data_mat, stim_mat_struc, stim_mat_simple, direc, frame_time)
%This function takes as inputs, the 3-D dff_data_mat that contains dF/F traces stored
%by frame_n, cell_n, trial_n and odor_n; and stim_mat which contains
%stimulus delivery information for each trial. The outputs are three 2D
%matrices. Resp_areas of size n_cell, n_trials contains the area under the pk response 
%of each cell on each trial sig_trace_mat of size n_cells, n_trials indicates which
%individual traces were significant responses, and sig_cell_mat of size 
%n_cells, n_odors which indicates which cells responded for more than half
%the presentations of any single odor.

n_frames = size(dff_data_mat, 1);
n_cells = size(dff_data_mat, 2);
n_trials = size(dff_data_mat, 3);
odor_list = unique(stim_mat_simple(:, 2));
n_odors = length(odor_list);
odor_dur_list = unique(stim_mat_simple(:, 3));
del = find(odor_dur_list == 0);
odor_dur_list(del) = [];
del = find(isnan(odor_dur_list) == 1);
odor_dur_list(del) = [];
n_odor_durs = length(odor_dur_list);

try
    stim_time = stim_mat_struc(1).stim_latency;              %stimulus onset time in s
catch
    stim_time = stim_mat_struc(1).stimLatency;              %stimulus onset time in s
end

stim_time = stim_time + 0.2;                                %added delay from valve opening to odor at pipe outlet can programmatically find this out - add feature later.
stim_frame = floor(stim_time./frame_time);                  %frame no at which odor reached fly
%pre_stim_frame = floor((stim_time - 4000)./frame_time);     %frame no 4 s prior to odor - to identify a 4s baseline
pre_stim_frame = 1;


resp_sizes = zeros(n_cells, n_trials) + nan;
resp_areaundercurves = zeros(n_cells, n_trials) + nan;
sig_trace_mat = zeros(n_cells, n_trials);
sig_cell_mat = zeros(n_cells, 12, n_odor_durs);

sig_trace_mat_old = zeros(n_cells, n_trials);
sig_cell_mat_old = zeros(n_cells, 12, n_odor_durs);


n_rep_gps = length(odor_list).*length(odor_dur_list);   %number of odor-dur pairs with sets of repeats
gp_list = [];

%builing list of odor, odor-dur pairs
for odor_n = 1:n_odors
    for odor_dur_n = 1:n_odor_durs
        gp_list = [gp_list; odor_list(odor_n), odor_dur_list(odor_dur_n)];
        
    end
end

for rep_gp = 1:n_rep_gps                %a rep_gp is a group of repeats of the same stimulus types
    curr_gp = gp_list(rep_gp, :);       %stimuli delivered in current group of repeats
    
    odor_ni = curr_gp(1);
    stim_duration = curr_gp(2);
    
    if isnan(stim_duration) == 1
        continue
    else
    end
    stim_end_fr = ceil(stim_frame + ((stim_duration.*1)./frame_time) );
    rep_tr_list = find(stim_mat_simple(:, 2) == odor_ni & stim_mat_simple(:, 3) == stim_duration);    %list of repeat tr numbers
    
    
    for cell_n = 1:n_cells
        curr_traces = squeeze(dff_data_mat(:, cell_n, rep_tr_list));             
        
        %moving window filtering each trace
        curr_traces_f = zeros(size(curr_traces, 1), size(curr_traces, 2)) + nan;
        for rep_n = 1:size(curr_traces, 2)
            curr_trace = curr_traces(:, rep_n);
            curr_traces_f(:, rep_n) = movmean(curr_trace, round(0.5./frame_time));        %filtering with a 0.5 s wide box-car
        end
        
       
        base_traces = curr_traces_f(pre_stim_frame:(stim_frame - 2), :);                 %pre-odor-stim baseline frames
        base_m = median(base_traces, 1, 'omitnan');
        base_s = std(base_traces, 1, 'omitnan');
        threshes = abs(base_m) + 3.*base_s;
        thresh_mean = mean(mean(base_traces, 2, 'omitnan'), 'omitnan') + 0.*std(mean(base_traces, 2, 'omitnan'), 'omitnan');       %significance criterion for resp in mean trace
        
        try
            resp_traces = curr_traces_f(stim_frame:round(stim_end_fr + 5./frame_time), :); %extending analysis window to 5s after odor off to capture off responses
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
                
            else
                
            end
            
            
            %in parallel checking if single trace peaks (not values at avg
            %trace pk) are significant as per the old criterion for 1s
            %stims
            if pk_resps_indiv(rep_n) > threshes(rep_n)
                sig_trace_mat_old(cell_n, rep_tr_list(rep_n)) = 1;
            else
                %keyboard
            end
                
        end
        
        %deciding if the current cell is a significant responder in the current rep-gp
        %if sum(sig_trace_mat(cell_n, rep_tr_list))./length(rep_tr_list) >= 0.5      %at least half the repeats should be sigonificant resps
            
        if sum(sig_trace_mat(cell_n, rep_tr_list)) > 2                          %at least 2 of the repeats shold be significant resps
                                                         
                dur_n = find(odor_dur_list == stim_duration);
                sig_cell_mat(cell_n, odor_ni, dur_n) = 1;
                h = 1;
%             else
%                 h = 0;
%             end
            
        else
            
            h = 0;
        end
        
        %repeating for old criterion
        if sum(sig_trace_mat_old(cell_n, rep_tr_list))./length(rep_tr_list) >= 0.5      %at least half the repeats should be sigonificant resps
            
            if sum(sig_trace_mat_old(cell_n, rep_tr_list)) > 2                          %at least 2 of the repeats shold be significant resps
                                                         
                dur_n = find(odor_dur_list == stim_duration);
                sig_cell_mat_old(cell_n, odor_ni, dur_n) = 1;
                h_old = 1;
                
            else
                h_old = 0;
            end
        else
            h_old = 0;
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
end

%sig_cell_1s_mat = squeeze(sig_cell_mat(:, :, 1));
