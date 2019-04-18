function [bad_trials] = magnet_plateau_detector(score_mat, beep_frame, puff_frame, frame_time)
%suntax: [bad_trials] = magnet_plateau_detector(score_mat)
%magnet_plateau_detector detects 'state-switches' which are present in a
%fair number of the magnet trace datasets. bad_trials is a vector with the 
%indices of all bad trials. A bad trial is a trial where there is a sudden 
%jump (state-switch) within the tone or trace periods. The hardware cause for these
%state-switches has subsequently been identified and fixed.

no_trials = size(score_mat, 2);
bad_trials = zeros(no_trials, 1);
for trial_no = 1:no_trials
    trace = score_mat(:, trial_no);
    
    %clipping away tone and trace periods for detection of state-switch (can't look in US period)
    detect_window = trace( (beep_frame - round(500./frame_time)):(beep_frame + (600./frame_time) ), 1);
    
    %down-sampling detect_window so as to detect peak in difference vector
    %as a signature of state-switch
    sub_s_factor = 150;     %factor by which vec is sub-sampled
    sub_s_vec = 1:sub_s_factor:(floor(length(detect_window)./sub_s_factor).*sub_s_factor);
    sub_s_vec = [sub_s_vec, length(detect_window)];
    trace_sub_s = detect_window(sub_s_vec);
    d_vec = diff(trace_sub_s);
    %saving detected state-switch as a bad trial
    if max(d_vec) > .4

          %USEFUL CODE FOR DEBUGGING/THRESHOLDING
%         figure(1)
%         subplot(2, 1, 1)
%         plot(trace_sub_s)
%         subplot(2, 1, 2)
%         plot(d_vec, '.')
%         keyboard

        d_pksi = find(abs(d_vec) > 0.2);
        if length(d_pksi) > 2       %legitimate peak, not DC offset
        else
            bad_trials(trial_no, 1) = 1;
            
        end
        
    else
    end
end


bad_trials = find(bad_trials == 1);
    