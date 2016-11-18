clear all
close all

stim_frame = 30;
stim_off_frame = 160;

tau1 = 10;
tau2 = 20;
ti = 1;
tf = 200;
t = linspace(ti,tf,200);
onset_vec1 = exp(-t/tau1).*3;
onset_vec2 = exp(-t/tau2).*3;

off_vec1 = circshift(onset_vec1', [stim_off_frame])';
off_vec2 = circshift(onset_vec2', [stim_off_frame])';

onset_vec1 = circshift(onset_vec1', [stim_frame])';
onset_vec2 = circshift(onset_vec2', [stim_frame])';

sus_vec1 = zeros(1, 200);
sus_vec1(stim_frame:stim_off_frame) = 3;

sus_vec2 = zeros(1, 200);
sus_vec2(stim_frame:stim_off_frame) = 2;

resp_mat = [onset_vec1; onset_vec2; off_vec1; off_vec2; sus_vec1; sus_vec2];
resp_mat_noise = resp_mat + randn(size(resp_mat, 1), size(resp_mat, 2));

resp_mat_smoothed = zeros(size(resp_mat, 1), size(resp_mat, 2));
lag = 10;
for trace_n = 1:size(resp_mat, 1);
    resp_mat_smoothed(trace_n, :) = tsmovavg_m(resp_mat_noise(trace_n, :)','s', lag, 1);
end

resp_mat_smoothed(:, 1:(lag-1)) = [];
 
dist_mat_clean =  corrcoef(resp_mat');
figure(1)
imagesc(dist_mat_clean)
title('clean traces')

dist_mat_noise = corrcoef(resp_mat_noise');
figure(3)
imagesc(dist_mat_noise);
title('noise added')

dist_mat_smoothed = corrcoef(resp_mat_smoothed');
figure(2)
imagesc(dist_mat_smoothed)
title('smoothed')
