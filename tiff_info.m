%by Mehrab N. Modi, 9/8/2014

function [no_frames, no_trials, frame_rate, zoom, d_path, f_name_out, tr_tg_no, ave_frame] = tiff_info(direc)

%direc = 'C:\Data\CSHL\20140908\no_stim_ave_im_after_zadjusted064.tif';

%identifying a filename in direc
dir_contents = dir([direc '*.tif']);
f_name = dir_contents(1, 1);
f_name = f_name.name;
f_name_out = f_name(1:(length(f_name) - 7) );
clear dir_contents

%identifying first trial's tag number
tr_tg_noi = findstr(f_name, '.tif');
tr_tg_no = str2double(f_name( (tr_tg_noi - 3):(tr_tg_noi - 1) ));



%reading in meta data
f_path = [direc f_name];
tif_tag = imfinfo([direc f_name]);
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

no_trials = dir([direc, '*.tif']);
no_trials = length(no_trials);

frame_ratei = findstr(img_desc, 'frameRate');
frame_rate = img_desc( (frame_ratei+10):(frame_ratei+16));
frame_rate = str2num(frame_rate);
clear frame_ratei


d_datei = findstr(direc, '201');
d_endi = findstr(direc, '\');
d_endi = d_endi(length(d_endi));

d_path = direc(d_datei:d_endi);
clear d_datei
clear d_endi


%reading in the first 30 frames of the first trial to obtain average image
for frame_no = 1:min([30, no_frames])
    frame = double(imread([direc f_name], frame_no));
    if frame_no == 1
        ave_frame = frame;
        continue
    else
    end
    ave_frame = ave_frame + frame;
end
ave_frame = ave_frame./30;




end






