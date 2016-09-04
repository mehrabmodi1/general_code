%by Mehrab N. Modi, 9/8/2014

direc = 'C:\Data\CSHL\20140908\no_stim_ave_im_after_zadjusted064.tif';

%reading in meta data
tif_tag = imfinfo(direc);
img_desc = tif_tag.ImageDescription;

zoomi = findstr(img_desc, 'zoomFactor');
zoom = img_desc( (zoomi+11):(zoomi+14));
zoom = str2num(zoom);
clear zoomi

no_framesi = findstr(img_desc, 'numberOfFrames');
no_frames = img_desc( (no_framesi+15):(no_framesi+18));
no_frames = str2num(no_frames);
%checking if no of frames is a 3 digit number
if isempty(no_frames) == 1
    i = 1;
    while isempty(no_frames) == 1
        i = 1+1;
        no_frames = img_desc( (no_framesi+15):(no_framesi+18 - i));
        no_frames = str2num(no_frames);
    end
else
end
clear no_framesi
clear i

no_trialsi = findstr(img_desc, 'numberOfRepeats');
no_trials = img_desc( (no_trialsi+16):(no_trialsi+19));
no_trials = str2num(no_trials);
%checking if no of frames is a 3 digit number
if isempty(no_trials) == 1
    i = 1;
    while isempty(no_trials) == 1
        i = 1+1;
        no_trials = img_desc( (no_trialsi+16):(no_trialsi+19 - i));
        no_trials = str2num(no_trials);
    end
else
end
clear no_trialsi
clear i


frame_ratei = findstr(img_desc, 'frameRate');
frame_rate = img_desc( (frame_ratei+10):(frame_ratei+16));
frame_rate = str2num(frame_rate);
clear frame_ratei






