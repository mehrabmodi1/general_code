%function [bad_frames] = bad_frame_finder(frame_mat, template_fr);

%This script thresholds frames just above background intensity and then
%compares the shape of the smoothed outline of the mushroom body in each
%frame with that of the averaged frame from early in the trial. It then identifies bad
%frames if too different from the reference frame.
%Mehrab Modi, 20141029


%=======================================
%loading in mock data for testing
frame_mat = load('C:\Data\CSHL\Analysed_data\sandbox\frame_mat2.mat');
frame_mat = frame_mat.saved_frame_mat;

template_fr = load('C:\Data\CSHL\Analysed_data\sandbox\template_fr2.txt');

%=======================================
%calculating outline for template frame
thresh = 150;           %obtained by manually examining frames

bin_template1 = binarise_im(template_fr, thresh);
bin_template2 = bwareaopen(bin_template1, 80);
bin_template = bwconvhull(bin_template2);

figure(1)
imagesc(bin_template)

thresh = 150;           %different thresh because these are single frames, while the template is averaged
%calculating outlines for each frame and comparing with template's outline
c_vec = zeros(1, size(frame_mat, 3)) + nan;
for frame_n = 1:size(frame_mat, 3)
    frame = frame_mat(:, :, frame_n);
    
    %replacing nans with 0's for the purpose of outline determination only (nans at edges due to rigid translation of frames)
    a = isnan(frame);
    a = find(a == 1);
    frame(a) = 0;
    
    bin_fr = binarise_im(frame, thresh);
    bin_fr = bwareaopen(bin_fr, 80);
    bin_fr = bwconvhull(bin_fr);
    
    c = (sum(sum(abs(bin_template - bin_fr))))./sum(sum(bin_template));     %the smaller this value, the better the match!
    c_vec(frame_n) = c;
    
    figure(2)
    imagesc(bin_fr)
    
    if frame_n == 44 | frame_n == 45 | frame_n == 67 | frame_n == 68
        keyboard
    else
    end
end



plot(c_vec, '.')
axis([0  100.0000    0    0.3200])

keyboard

bad_frames = [];

%end