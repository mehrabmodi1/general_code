clear all
close all


path = 'C:\Data\Data\Raw_data\20180508\';
cd(path);
tif_list = dir('*.tif');

n_tifs = size(tif_list, 1);

scale = 704;        %in pixels/mm

for tif_n = 1:n_tifs

    info = imfinfo([path, tif_list(tif_n).name]);
    n_frames = size(info, 1);
    
    %skipping scale image
    if n_frames == 1
        continue
    else
    end


    %loop to read in all frames
    dia_vec = [];
    for frame_n = 1:n_frames
        curr_frame = imread([path, tif_list(tif_n).name], frame_n);
        curr_frame_cropped = curr_frame(415:531, 549:871);

        %binarising
        curr_frame_cropped(curr_frame_cropped<5) = 0;
        curr_frame_cropped(curr_frame_cropped>=5) = 1;

        bnd_boxes = regionprops(curr_frame_cropped, 'BoundingBox');
        dia = bnd_boxes(1).BoundingBox(3);
        dia_vec = [dia_vec; dia];

    end
    
    dia_vec = dia_vec./scale;
    rand_vec = rand(n_frames, 1) + (tif_n - 1).*1.5;
    sd = std(dia_vec)
    
    figure(1)
    plot(rand_vec, dia_vec, 'o')
    ylabel('ball diameter (mm)')

    axis_vals = axis;
    axis_vals = [-1, (tif_n.*1.5), 0, axis_vals(4).*1.2];
    axis(axis_vals);
    hold on
    drawnow
    
    

end
