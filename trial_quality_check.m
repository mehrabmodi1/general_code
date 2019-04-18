%function [data_mat_corrected] = trial_quality_check(dataset, trial_n)
%This function identifies bad frames as those that 

clear all
close all
dataset = load('C:\Data\CSHL\20150428\expt\dopamine_expt.mat');
dataset = dataset.data;
trial_n = 1;

within_trial_thresh = 0.49;

for trial_n = 1:size(dataset, 2)

    raw_stack = dataset(trial_n).imageStack;
    n_frames = size(raw_stack, 3);
    stim_time = dataset(trial_n).stim.stimLatency;
    frame_time = dataset(trial_n).info.framePeriod;
    stim_frame = floor(stim_time./frame_time);
    stim_frames = [(stim_frame + round(1.15./frame_time)), (stim_frame + round(4./frame_time))];   %ignoring 1s after odor arrival bec fly usually moves then              
    n_s_frames = stim_frames(1, 2) - stim_frames(1, 1);

    %calculating difference frames and then amplifying differences with
    %imopen
    counter = 0;
    counteri = 0;
    corr_mat = zeros(n_s_frames, n_s_frames);
    for frame_n = stim_frames(1, 1):stim_frames(1, 2)
        counter = counter + 1;
        counteri = counter;

        for frame_ni = frame_n:stim_frames(1, 2)
            counteri = counteri + 1;
            frame1 = raw_stack(:, :, frame_n);
            frame2 = raw_stack(:, :, frame_ni);
            diff_frame = abs(frame2 - frame1);
            diff_frame1 = imtophat(diff_frame, strel('disk', 5));
            diff_score = sum(sum(diff_frame));
            
            figure(1)
            subplot(2, 2, 1)
            imagesc(frame1);
            title(['dff_score = ' int2str(diff_score)])
            colormap('gray')
            subplot(2, 2, 2)
            imagesc(frame2);
            subplot(2, 2, 3)
            imagesc(diff_frame);
            subplot(2, 2, 4)
            imagesc(diff_frame1);
            keyboard

        end


    end

   
    
    keyboard
end

%end