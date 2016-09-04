function [] = trial_movt_chk(direc, trial)
%This function compares the pre-stim averaged frames of trial with the
%first trial in direc. It also displays the two averaged images. The
%comparison is made by calculating the correlation coefficients for the
%5-frame and all available frame averages.

[n_frames, n_trials, frame_rate, zoom, d_path, f_name, tr_tg_num, ave_frame] = tiff_info(direc);
[isi, odor_t_list, odorNames, n_reps, stim_time] = load_stim_param(direc);

frame_time = 1./frame_rate;                                             %in s
stim_frame = floor(stim_time(1, 1)./frame_time);                              %frame number at which stim was delivered

%loop to load frames for the first trial and the test trial
for trial_n = 1:2
    
    %loop to load frames
    for frame_n = 1:stim_frame
        
    end
    
end



end