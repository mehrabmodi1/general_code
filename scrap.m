clear all
close all

im = zeros(100, 100);
im(30:70, 40:60) = 2;
im2 = zeros(100, 100);
im2(30:70, 40:60) = 2;
im2(15:60, 10:35) = 2;

r_shifts = round(rand(2, 50).*25);
r_stack = zeros(100, 100, 100);
for frame_n = 1:50
    curr_im = circshift(im, r_shifts(1, frame_n), 1);
    curr_im = circshift(curr_im, r_shifts(2, frame_n), 2);
    r_stack(:, :, frame_n) = curr_im + rand(100, 100);
    curr_im = circshift(im2, r_shifts(1, frame_n), 1);
    curr_im = circshift(curr_im, r_shifts(2, frame_n), 2);
    r_stack(:, :, (frame_n + 50)) = curr_im + rand(100, 100);
end


test_frame = r_stack(:, :, 55).*0.5;
c_vec = [];
c_vec2 = [];
for frame_n = 1:100
    c = mat_corrcoef(test_frame, r_stack(:, :, frame_n));
    c_vec = [c_vec; c];
    
    test_frame = (test_frame - mean(mean(test_frame)));
    %test_frame = test_frame./max(max(test_frame));
    ref_frame = r_stack(:, :, frame_n) - mean(mean(r_stack(:, :, frame_n)));
    %ref_frame = ref_frame./max(max(ref_frame));
    c2_mat = xcorr2(test_frame, ref_frame);
    c_vec2 = [c_vec2; max(max(c2_mat))];
    
end

figure(1)
plot(c_vec)
hold on
%plot(c_vec2, 'r')