function [diff_fr, sig_resp_frame, sig_resp_frame_n, smple_fr] = highlight_resp_pix(fig_n, stack, stim_frs, n_channels, frame_time, suppress_plot)
%syntax: [diff_fr, sig_resp_frame, sig_resp_frame_n, smple_fr] = highlight_resp_pix(fig_n, stack, stim_frs, frame_time, suppress_plot)
%This function identifies pixels whose dF/F response is in the top
%cutoff_pc percent of all pixels in the mean response frame. It then plots
%the mean response frame with the responsive pixels highlighted.


%identifying responsive pixels
resp_fr_final = zeros(size(stack, 1), size(stack, 2));

if n_channels == 2
    stack = stack(:, :, 1:2:size(stack, 3));
else
end

%for pulse_n = 1:size(stim_frs, 1)
for pulse_n = 1:size(stim_frs, 1)
    stim_fr = stim_frs(pulse_n, 1);
    stim_end_fr = stim_frs(pulse_n, 2);
    baseline_frs = stack(:, :, (stim_fr - round(4./frame_time) ):(stim_fr - 1) );
    baseline_frs = double(baseline_frs);
    if pulse_n == 1
        baseline_fr = mean(baseline_frs, 3);
        baseline_sds = std(baseline_frs, 0, 3);
        bk_pixi = (baseline_fr < 7);
        baseline_fr(bk_pixi) = 0;
        baseline_sds(bk_pixi) = 0;
    else
    end
    
    resp_fr_curr = mean(stack(:, :, (stim_fr + 1):(stim_end_fr + round(2./frame_time)) ), 3);
    resp_fr_final = resp_fr_final + resp_fr_curr;
end
resp_fr = resp_fr_final./size(stim_frs, 1);       %averaged response frame across all stimulus pulses
resp_fr(bk_pixi) = 0;
resp_fr = double(resp_fr);
    
diff_fr = (resp_fr - baseline_fr);
diff_fr = diff_fr - 2.33.*baseline_sds;               %identifying pixels with a dF > 2.3 SDs of baseline period fluctuations ie significantly responsive
diff_fr = diff_fr./baseline_fr;                      %dividing by F to get a dF/F image
diff_fr(diff_fr < 0) = 0;                            %getting rid of negs
nansi = isnan(diff_fr);
diff_fr(nansi) = 0;                                  %getting rid of nans
diff_fr(diff_fr > 10) = 10;                          %forcing crazy dF/F values to sane ones
diff_fr = imgaussfilt(diff_fr, 2);

bin_pix = diff_fr > 0;
bin_pix = bwareaopen(bin_pix, 100);                  %getting rid of small objects
diff_fr = diff_fr.*bin_pix;


%getting rid of background pixels
bin_pix(baseline_fr < 7) = 0;
diff_fr = diff_fr.*bin_pix;
smple_fr(baseline_fr > 7) = 1;

sig_resp_frame = diff_fr.*bin_pix;

% %test plots
% figure(1)
% imagesc(baseline_fr_n)
% figure(2)
% imagesc(diff_fr, [-3, 4])

%plotting
cust_map = colormap([gray(64); autumn(64)]);
resp_fr_n = resp_fr./max(max(resp_fr));
diff_fr_n = diff_fr./max(max(diff_fr));
%making sure to split the colormap by ofsetting diff_fr values by 50%
resp_fr(resp_fr_n < 0) = 0;     %ranges from 0 to 1
diff_fr_n = diff_fr_n + 1;      %ranges from 1 to 2

if suppress_plot == 0
    figure(fig_n)
    colormap(cust_map)
    base_im = imagesc(resp_fr_n);
    hold on
    s = imagesc(diff_fr_n);
    alphamap = ((diff_fr_n - 1).*bin_pix).*1;
    alpha(s, alphamap)
    hold off
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
else
end

sig_resp_frame_n = (diff_fr_n - 1).*bin_pix;

