function n_frames = parse_nframes(im_path)
%This function takes the path string for an SI-saved multi-page tiff file.
%It parses the ImageDescription meta-data to find out the number of frames
%in the file.
%Mehrab Modi, 20160511

info = imfinfo(im_path);
info = info(1).ImageDescription;
del = strfind(info, 'numberOfFrames');
n_frames_str = [];
for chari = 1:6
    curr_char = info((del+14) + chari);
    if isempty(str2num(curr_char)) == 1
        break
    else
        n_frames_str = [n_frames_str, curr_char];
    end
end
n_frames = str2num(n_frames_str);

end