function [raw_data_mat no_frames CS_onset_frame US_onset_frame] = rand_raw_data_clipper(raw_data_mat, set_type,...
    stimOI, rand_times_list, pre_clip, post_clip, frame_time, CS_onset_frame, US_onset_frame, CS_US_delay)
%syntax: [raw_data_mat no_frames CS_onset_frame US_onset_frame] =
%rand_raw_data_clipper(raw_data_mat, pre_clip, post_clip, rand_times_list, set_type)
%This function clips away data farther away from the tone/puff periods than specified by pre_clip
%and post_clip, taking into account whether or not the dataset is a trace, a control or a rand
%dataset. The function pads with Nans if the tone/puff came too close to
%the beginning or end of the trial.

if nargin == 1
    set_type = 0;
    stimOI = 1;
    rand_times_list = [];
    pre_clip = 2000;
    post_clip = 2500;
    frame_time = 80;    %ms
elseif nargin == 2
    stimOI = 1;
    rand_times_list = [];
    pre_clip = 2000;
    post_clip = 2500;
    frame_time = 80;    %ms
elseif nargin == 3
    rand_times_list = [];
    pre_clip = 2000;
    post_clip = 2500;
    frame_time = 80;    %ms
elseif nargin == 4
    pre_clip = 2000;
    post_clip = 2500;
    frame_time = 80;    %ms
elseif nargin == 5
    post_clip = 2500;
    frame_time = 80;    %ms
elseif nargin == 6
    frame_time = 80;    %ms
end

no_frames = size(raw_data_mat, 1);
no_cells = size(raw_data_mat, 2);
no_trials = size(raw_data_mat, 3);

pre_tone_frame_no = floor(pre_clip./frame_time);
post_tone_frame_no = floor(post_clip./frame_time);

%CLIPPING data near point of interest (tone or puff) 
no_framesa = pre_tone_frame_no + post_tone_frame_no + 1;
raw_data_mati = zeros(no_framesa, no_cells, no_trials);
%clear no_framesa
clear tone_framea
if set_type == 0                     %case when dataset is trace or control_no-puff
    raw_data_mati = raw_data_mat( (CS_onset_frame - pre_tone_frame_no):(CS_onset_frame + post_tone_frame_no), :, :);
    raw_data_mat = raw_data_mati;
    clear raw_data_mati
elseif set_type == 1 && stimOI == 1  %case when dataset is rand and stimulus of interest is tone
    for trial_no = 1:no_trials
        tone_frame = round(rand_times_list(1, trial_no)./frame_time);

        %checking if rand tone came too close to time = 0 or end, padding to align with other trials
        if floor(tone_frame - pre_tone_frame_no)>0 && floor(tone_frame + post_tone_frame_no)<no_frames     %case where tone came far from 0 and end   
            raw_data_mati(:, :, trial_no) = raw_data_mat((tone_frame - pre_tone_frame_no):(tone_frame + post_tone_frame_no), :, trial_no);
        elseif floor(tone_frame - pre_tone_frame_no)<0                                                     %case where tone came close to 0
            pad = zeros(abs(floor(tone_frame - pre_tone_frame_no)) + 1, no_cells) + nan;                   %0's replaced with nans so that they dont contribute to PSTHs (calculated using nanmean)
            raw_data_mati(:, :, trial_no) = [pad; squeeze(raw_data_mat( 1:floor(tone_frame + post_tone_frame_no), :, trial_no))];
        elseif floor(tone_frame + post_tone_frame_no)>no_frames                                             %case where tone came close to end
            pad = zeros((floor(tone_frame + post_tone_frame_no) - no_frames), no_cells) + nan;              %0's replaced with nans so that they dont contribute to PSTHs
            raw_data_mati(:, :, trial_no) = [squeeze(raw_data_mat(floor(tone_frame - pre_tone_frame_no):end, :, trial_no)); pad];
        end
        clear trial_no
    end

    raw_data_mat = raw_data_mati;
    clear raw_data_mati
    clear tone_frame
elseif set_type == 1 && stimOI == 2  %case when dataset is rand and stimulus of interest is puff
    raw_data_mat = raw_data_mat( round(US_onset_frame - pre_tone_frame_no):round(US_onset_frame + post_tone_frame_no), :, :);
end
%re-assigning new frame numbers for imp events post-clipping
no_frames = size(raw_data_mat, 1);
CS_onset_frame = pre_tone_frame_no;
US_onset_frame = CS_onset_frame + round(CS_US_delay./frame_time);
