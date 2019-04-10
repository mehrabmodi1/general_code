function [r_vecs_saved, h_vecs_saved] = cell_classifier(cell_data, frame_time)
all_traces = cell_data.traces;
n_odors = size(all_traces, 2);
n_odor_durs = size(all_traces, 3);
stim_frs = [cell_data.stim_start_frs, cell_data.stim_end_frs];

sig_thresh = 0.01;

manual_assessment = 1;

r_vecs_saved = zeros(n_odors, 6) + nan;
h_vecs_saved = zeros(n_odors, 6) + nan;

%analysis over entire stimulus window
for odor_n = 1:n_odors
    h_vec_mat = [];
    r_vec_mat = [];
    follow_up = 0;
    for odor_dur_n = 1:n_odor_durs
        curr_trace_m = nanmean(squeeze(cell_data.traces(:, odor_n, odor_dur_n, :) ), 2);
        
        %for lost trials
        if sum(isnan(curr_trace_m)) == length(curr_trace_m)
            disp(['All trials lost for this cell, for odor_dur', int2str(odor_dur_n) ' odor ' int2str(odor_n)]);
            h_vec = zeros(1, 5) + nan;
            result_vec = h_vec;
            %keeping track of h_vec across odor durations
            h_vec_mat = [h_vec_mat; h_vec];
            r_vec_mat = [r_vec_mat; result_vec];
            %keyboard
            continue
        else
        end
        
        n_frames = length(curr_trace_m);
        start_fr = stim_frs(odor_dur_n, 1);
        end_fr = stim_frs(odor_dur_n, 2);
        odor_dur_ni = (end_fr - start_fr).*frame_time;      %odor duration in s
        
        stim_win = [start_fr, end_fr];
        stim_trace = curr_trace_m(stim_win(1, 1):stim_win(1, 2));
        
                
        %smoothing win trace to detect number of threshold crossings
        filt_width = 5;
        filter_ = zeros(1, filt_width) + 1./filt_width;
        s_stim_trace = filter(filter_, 1, stim_trace);
                
        %setting stim_win cross threshold, counting threshold crossings
        s_thresh = max(s_stim_trace).*0.25;
        del = s_stim_trace < s_thresh;      
        s_stim_trace_th = ones(1, length(s_stim_trace));
        s_stim_trace_th(del) = 0;           %thresholded, smoothed curve
        s_stim_trace_th = bwareaopen(s_stim_trace_th, 10);   %getting rid of noisy, rapid fluctuations near threshold
        
        
        figure(5)
        imagesc(s_stim_trace_th)
        figure(6)
        plot(s_stim_trace)
        label_mat = bwlabel(s_stim_trace_th);
        
        if max(label_mat) > 1
            for pk_n = 1:max(label_mat)
                curr_pki = find(label_mat == pk_n);
                trace_seg = s_stim_trace(curr_pki);                      %calculating areas above threshold of each detected peak
                area_vec(pk_n) = sum(trace_seg) - (s_thresh.*length(trace_seg) );
            end
            
            max_pk = max(area_vec);
            small_pks = find(area_vec < (max_pk.*0.7) );
            area_vec(small_pks) = 0;
            big_pks = find(area_vec > 0);
            
            if length(big_pks) > 1
                multi_pk = 1;       %This cell is a complex response cell with multiple peaks in the odor window
            else
                multi_pk = 0;
            end
            r_multi_pk = mean(area_vec(big_pks));
            
        else
            multi_pk = 0;
            r_multi_pk = 0;
        end
        
        
        %--------------------------------------------------
        
        
        %splitting response curve into smaller windows to do a finer
        %analysis
        on_width = odor_dur_ni.*0.3;       %width of on-cell analysis window, beginning from stim on, in s
        on_width_s = 1.5;                  %width of short on-cell analysis window, beginning from stim on, in s 
        ramp_width_s = 2.5;                  %width of short ramp region window, only for short on window comparison     
        sus_width = odor_dur_ni.*0.3;      %width of sus-cell analysis window, going backwards from stim off, in s
        off_width = odor_dur_ni.*0.3;      %width of off-cell analysis window, starting from stim off, in s
        off_width_s = 1.5;                    %width of short off-cell analysis window, beginning from stim on, in s
        ramp_width = odor_dur_ni.*0.4;     %width of ramp-cell analysis window, going backwards from beginning of sustained window, in s
        base_width = 2;     %width of baseline analysis window, going backwards from stim on, in s
                
        %skipping step if stim duration was longer than acquisition time
        if end_fr > n_frames
            
            h_vec = zeros(1, 5) + nan;
            result_vec = zeros(1, 5) + nan;
            %keeping track of h_vec across odor durations
            h_vec_mat = [h_vec_mat; h_vec];
            r_vec_mat = [r_vec_mat; result_vec];
            continue
        else
        end
        
        
        
        on_win = [start_fr, (start_fr + round(on_width./frame_time) )];
        base_win = [(start_fr - round(base_width./frame_time)), (start_fr - 1)];
        
        %skipping current stim duration if its too short for non-overlapping
        %on-ramp-sustained-off analysis windows
        if  (on_width + sus_width + ramp_width) < 1.2       %three 400 ms Ca decays
            h_vec = zeros(1, 5) + nan;
            result_vec = h_vec;
            
            on_win = [start_fr, (start_fr + round(1.5./frame_time) )];
            base_win = [(start_fr - round(base_width./frame_time)), (start_fr - 1)];
                        
            %calculating on response size
            on_trace = curr_trace_m(on_win(1, 1):on_win(1, 2));
            base_trace = curr_trace_m(base_win(1, 1):base_win(1, 2));
            on_diff = nanmean(on_trace) - nanmean(base_trace);
            result_vec(1, 1) = on_diff;
            
            %statistical testing for significance of response
            [h_on, p_on] = ttest2(base_trace, on_trace, 'Alpha', sig_thresh);
            h_vec(1, 1) = h_on;
            
            %keeping track of h_vec across odor durations
            h_vec_mat = [h_vec_mat; h_vec];
            r_vec_mat = [r_vec_mat; result_vec];

            continue
        else
        end
       
        
        on_win_s = [start_fr, (start_fr + round(on_width_s./frame_time) )];
        ramp_win_s = [(on_win_s(1, 2) + 1), ( on_win_s(1, 2) + round(ramp_width_s./frame_time) )];
        sus_win = [(end_fr - round(sus_width./frame_time) ), end_fr];
        off_win = [(end_fr + 1), (end_fr + round(off_width./frame_time) )];
        off_win_s = [(end_fr + 1), (end_fr + round(off_width_s./frame_time) )];
        ramp_win = [(sus_win(1, 1) - round(ramp_width./frame_time) ), (sus_win(1, 1) - 1)];
        
        
        
        if off_win(1, 2) > n_frames
            off_win(1, 2) = n_frames;
        else
        end
        
        on_trace = curr_trace_m(on_win(1, 1):on_win(1, 2));
        on_trace_s = curr_trace_m(on_win_s(1, 1):on_win_s(1, 2));
        ramp_trace_s = curr_trace_m(ramp_win_s(1, 1):ramp_win_s(1, 2));
        off_trace = curr_trace_m(off_win(1, 1):off_win(1, 2));
        off_trace_s = curr_trace_m(off_win_s(1, 1):off_win_s(1, 2));
        sus_trace = curr_trace_m(sus_win(1, 1):sus_win(1, 2));
        ramp_trace = curr_trace_m(ramp_win(1, 1):ramp_win(1, 2));
        base_trace = curr_trace_m(base_win(1, 1):base_win(1, 2));
        
        
        %calculating differences between means 
        on_diff = nanmean(on_trace) - nanmean(base_trace);
        on_ramp_diff = nanmean(on_trace) - nanmean(ramp_trace);
        off_diff = nanmean(off_trace) - nanmean(sus_trace);
        sus_diff = nanmean(sus_trace) - nanmean(base_trace);
        ramp_diff = nanmean(sus_trace) - nanmean(ramp_trace);
        
        
        %identifying windows for this trace that are significantly higher
        %than the baseline window
        [h_on, p_on] = ttest2(base_trace, on_trace, 'Alpha', sig_thresh);
        [h_on2, p_on2] = ttest2(ramp_trace, on_trace, 'Alpha', sig_thresh);
        [h_on_s, p_on_s] = ttest2(base_trace, on_trace_s, 'Alpha', sig_thresh);            %for shorter on window
        [h_on2_s, p_on2_s] = ttest2(ramp_trace, on_trace_s, 'Alpha', sig_thresh);          %for shorter on window 
        
        [h_off, p_off] = ttest2(sus_trace, off_trace, 'Alpha', sig_thresh);
        [h_off_s, p_off_s] = ttest2(sus_trace, off_trace_s, 'Alpha', sig_thresh);

        [h_sus, p_sus] = ttest2(base_trace, sus_trace, 'Alpha', sig_thresh);
        [h_ramp, p_ramp] = ttest2(sus_trace, ramp_trace, 'Alpha', sig_thresh);
        
        %rules for determining if a difference is significant:
        %on responses should be higher than baseline window AND ramp window
        %ramp responses should be higher than baseline and lesser than sus;
        %but also, sus should be higher than baseline
        h_vec = [max([h_on.*h_on2, h_on_s.*h_on2_s]), h_ramp.*h_sus, h_sus, max([h_off, h_off_s])];
               
        %logic conditions on signs of differences to throw out spurious
        %catches
        %1. Off responses have to rise above the sustained level
        if sign(off_diff) == -1
            h_vec(1, 4) = 0;
        else
        end
        
        %2. Responses in the ramp window have to be lower than the
        %sustained level
        if sign(ramp_diff) == -1
            h_vec(1, 2) = 0;
        else
        end
            
        %3. Sustained response has to be higher than baseline
        if sign(sus_diff) == -1
            h_vec(1, 3) = 0;
        else
        end
        
        %4. On response has to be higher than baseline
        if sign(on_diff) == -1
            h_vec(1, 1) = 0;
        else
        end
        
        %5. On response has to be higher than sustained response
        if sign(on_ramp_diff) == -1
            h_vec(1, 1) = 0;
        else
        end
        
         
        if multi_pk == 1
            h_vec = [zeros(1, 4), 1];
        elseif multi_pk == 0
            h_vec = [h_vec, 0];
        end
        result_vec = [on_diff, ramp_diff, sus_diff, off_diff, r_multi_pk];         %vector of effect sizes, non zero only if significantly different
        
                       
        if manual_assessment == 1
            figure(odor_dur_n)
            %plotting for manual inspection of classification
            patch_wins = [on_win; ramp_win; sus_win; off_win];
            patch_colors = repmat([0.75, 0.5, 0.35], 4, 1);
            patch_colors = patch_colors .* repmat(h_vec(1:4)', 1, 3);
            plot(curr_trace_m, 'LineWidth', 2)
            title(['Multi-Pk' int2str(h_vec(5))]);
            xlabel('frame_n')
            ylabel('dF/F')
            add_stim_shading(odor_dur_n, patch_wins, 0.25, patch_colors)
            set_xlabels_time(odor_dur_n, frame_time, 1)
                
%             if odor_dur_n == 3
                %del = input('Press Enter')
%             else
%             end
            
        else
        end
        
        %keeping track of h_vec across odor durations
        h_vec_mat = [h_vec_mat; h_vec];
        r_vec_mat = [r_vec_mat; result_vec];
        
    end
    
    %reconciling differences in classification across multiple durations,
    %if any by asking user to decide
    if size(h_vec_mat, 1) > 1
        m_diff = 0;
        for n_h_vec = 2:size(h_vec_mat, 1)
            c_diff = sum(abs(nanmean(h_vec_mat, 1) - h_vec_mat(n_h_vec, :)));
            m_diff = max([m_diff, c_diff]);
        end
        
        long_h_vec = h_vec_mat(n_odor_durs, :);
        if long_h_vec(5) == 0
            h_vec = long_h_vec;
            m_diff = 0;
        else
            h_vec = h_vec_mat((n_odor_durs - 1), :);
            m_diff = 0;
        end
        
        
        %creating dialog box for user to intervene if classifier results were different
        %for two durations of the same odor for the same cell
        if m_diff > 0
            autoArrangeFigures(0, 0);
                        
            qstring = 'Which fig has correct windows identified as significant?';
            close(figure(1))
            for fig_n = 1:size(h_vec_mat, 1)
                btn_options{1, fig_n} = ['fig' int2str(fig_n+1)];
            end
            
            fig_n = bttnChoiseDialog(btn_options, 'Input Needed', 1, qstring);
            h_vec = h_vec_mat(fig_n, :);
            
        else
            h_vec = h_vec_mat(n_odor_durs, :);
        end
    
    else
    end
    
    try
        r_vec = r_vec_mat(size(r_vec_mat, 1), :);       %taking response sizes for the longest duration stimulus
        if sum(isnan(r_vec)) == length(r_vec)           %if longest duration stimulus trials were missing
            r_vec = nanmean(r_vec_mat(2:size(r_vec_mat, 1), :), 1);
        else
        end        
        
    catch
        keyboard
    end
    
    %saving shortest stim duration response size as well
    r_vec = [r_vec, r_vec_mat(1, 1)];       
    h_vec = [h_vec, h_vec_mat(1, 1)];
    
    r_vecs_saved(odor_n, :) = r_vec;
    h_vecs_saved(odor_n, :) = h_vec;
        
end

close(figure(3))
end