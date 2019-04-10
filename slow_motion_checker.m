function [bad_trials] = slow_motion_checker(direc, save_direc)
%this code is for data quality control. It calculates a 5-trial average at the
%beginning of each trial and compares this with the first trial by
%calculating a correlation coefficient. A line is fitted to all trial
%corrcoefs, a threshold line is calculated based on manually compared
%averaged images and bad trials identified based on the threshold.
%Mehrab Modi, 20141003


ave_c_vec = [];

%House-keeping
%-----------------------

direc = 'C:\Data\CSHL\20141028\ag_odor_stim\';                        %directory with raw data in it. Don't forget the last '\'!
n_aved = 5;                                                             %number of frames averaged to generate frames for comparison within trial 1 and between tr 1 and other trials


[n_frames, n_trials, frame_rate, zoom, d_path, f_name, tr_tg_num, ave_frame] = tiff_info(direc);
[isi, odor_t_list, odorNames, n_reps, stim_time] = load_stim_param(direc);

frame_time = 1./frame_rate;                                             %in s
stim_frame = floor(stim_time(1, 1)./frame_time);                              %frame number at which stim was delivered


ave_mat = zeros(size(ave_frame, 1), size(ave_frame, 2), floor(stim_frame./n_aved));

for trial_n = 1:n_trials
    %building filename string
    trial_no_f = trial_n + tr_tg_num - 1;           %adding the trial tag number of the first trial of this set to get the right filename
    trial_no_f = sprintf('%03.0f', trial_no_f);
    file_path = [direc f_name trial_no_f, '.tif'];  %full file path for currently analysed trial

    %building distribution of expected correlation coefficients for no
    %motion (comparing averaged frames within first trial's pre-stim period)
    
    %reading in 3 frame segments to calculate a matrix of averaged images
    frame_counter = 1;
    for frame_n = 1:n_aved:(stim_frame - rem(stim_frame, n_aved))

        for sub_frame_n = 0:(n_aved - 1)
            frame = double(imread(file_path, (frame_n + sub_frame_n) ) );

            if sub_frame_n == 0
                ave_frame = frame;
            else
            end

            if sub_frame_n > 0
                ave_frame = ave_frame + frame;
            else
            end
        end

        ave_frame = ave_frame./n_aved;
        ave_mat(:, :, frame_counter) = ave_frame;           % matrix storing averaged frame sets obtained after averaging 5 frames
        frame_counter = frame_counter + 1;
        
    end

    if trial_n == 1
        %calculating distribution of corrcoefs for random pairs of averaged
        %frames from within trial 1.
        ave_mat1 = ave_mat;
        template_frame = ave_mat(:, :, 1);                  %template frame from trial 1 compared with all other trials
        n_ave_frames = size(ave_mat, 3);
                       
    else
    end
    
    %comparing each frame in ave_mat with the template frame from trial 1
    curr_c_vec = zeros(size(ave_mat, 3), 1) + nan;
    for ave_f_n = 1:size(ave_mat, 3)
        curr_ave_f = ave_mat(:, :, ave_f_n);
        curr_c_vec(ave_f_n, 1) = mat_corrcoef(curr_ave_f, template_frame);
        
    end
    
    ave_c_vec = [ave_c_vec; mean(curr_c_vec)];
    
    
    
    %fitting a line to the corrcoef decay curve
    x = (1:1:length(ave_c_vec))';                       %vector for fitting and line calculation
    a = polyfit(x, ave_c_vec, 1);
    line = polyval(a, x);
    
end

figure(75)
plot(ave_c_vec, '.')
hold on
plot(line, 'r')
xlabel('trial number')
ylabel('corrcoef with trial 1')
title('Slow motion checker')

set(gcf, 'Color', 'w')

%thresholding at a corr-coef value below the fitted line to identify bad
%trials
thresh = line(1, 1).*0.0178;         %this multiplier was estimated manually, using data from C:\Data\CSHL\20140908\odor_stim_set2\
thresh_line = line - thresh;

%comparing theshold line with actually observed corrcoefs
com_mat = [ave_c_vec, thresh_line]';
[del, bad_trials] = min(com_mat);
bad_trials = find(bad_trials == 1);

clear del

end