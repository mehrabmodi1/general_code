function [diff_fr, sig_resp_frame, sig_resp_frame_n, baseline_fr_orig, stack] = highlight_resp_pix(fig_n, stack, stim_frs, n_channels, frame_time, pix_per_mic, suppress_plot)
%syntax: [diff_fr, sig_resp_frame, sig_resp_frame_n, smple_fr] = highlight_resp_pix(fig_n, stack, stim_frs, n_channels, frame_time, suppress_plot)
%This function identifies pixels whose dF/F response is in the top
%cutoff_pc percent of all pixels in the mean response frame. It then plots
%the mean response frame with the responsive pixels highlighted.

%converting stack to a double-precision variable
if isa(stack, 'double') == 0
    stack = double(stack);
else
end


if n_channels == 2
    stack = stack(:, :, 1:2:size(stack, 3));
else
end

%computing mean and sd of baseline frames
try
    baseline_frs = stack(:, :, (stim_frs(1,1) - round(4./frame_time) ):(stim_frs(1, 1) - 1) );
catch
    keyboard
end

baseline_fr = mean(baseline_frs, 3);

%computing background value of image (PMT offset)
bk_offset = mean(mean(baseline_fr(1:10, 1:10), 'omitnan'), 'omitnan');
bk_pixi = (baseline_fr < bk_offset);
stack = stack - bk_offset;

baseline_fr = mean(baseline_frs, 3);
baseline_sds = std(baseline_frs, [], 3);

%spatially averaging baseline sds 
se = strel('disk', 6);
baseline_sds = imopen(baseline_sds, se);

%computing response frame averaged across stim pulses
resp_fr_final = zeros(size(stack, 1), size(stack, 2));
for pulse_n = 1:size(stim_frs, 1)
    stim_fr = stim_frs(pulse_n, 1);
    stim_end_fr = stim_frs(pulse_n, 2);
    
    resp_fr_curr = mean(stack(:, :, (stim_fr + 1):(stim_end_fr + round(2./frame_time)) ), 3, 'omitnan');
    resp_fr_final = resp_fr_final + resp_fr_curr;
end

resp_fr = resp_fr_final./size(stim_frs, 1);       %averaged response frame across all stimulus pulses
    
diff_fr = (resp_fr - baseline_fr);

diff_fr(diff_fr < 2.*baseline_sds) = 0;               %identifying pixels with a dF > 2.3 SDs of baseline period fluctuations ie significantly responsive
diff_fr = diff_fr./baseline_fr;                      %dividing by F to get a dF/F image
diff_fr(diff_fr < 0) = 0;                            %getting rid of negs
nansi = isnan(diff_fr);
diff_fr(nansi) = 0;                                  %getting rid of nans
diff_fr(diff_fr > 10) = 10;                          %forcing crazy dF/F values to sane ones
diff_fr = imgaussfilt(diff_fr, 2);

%getting rid of small objects
bin_pix = diff_fr > 0;
bin_pix = bwareaopen(bin_pix, 40);                  
diff_fr = diff_fr.*bin_pix;

%getting rid of background pixels
% bin_pix(bk_pixi) = 0;
% diff_fr = diff_fr.*bin_pix;
sig_resp_frame = diff_fr.*bin_pix;

% %test plots
% figure(1)
% imagesc(baseline_fr_n)
% figure(2)
% imagesc(diff_fr, [-3, 4])

%plotting
cust_map = colormap([gray(512); autumn(512)]);
baseline_fr_orig = baseline_fr;
baseline_fr_plt = mean(stack, 3, 'omitnan');
%baseline_fr_plt = real(sqrt(baseline_fr_plt));
resp_fr_n = baseline_fr_plt./(quantile(reshape(baseline_fr_plt, 1, []), .99));
resp_fr_n(resp_fr_n >= 1) = 0.999;

pk_dff = max(max(diff_fr));
diff_fr_n = diff_fr./pk_dff;
if sum(sum(diff_fr_n)) > 0
else
    disp('Warning: No response pixels found. Highlighting baseline image instead.')
end

%making sure to split the colormap by offsetting diff_fr values by 50%
resp_fr_n(resp_fr_n < 0) = 0;     %ranges from 0 to 1
diff_fr_n = diff_fr_n + 1;      %ranges from 1 to 2
resp_fr_n = add_scale_bar(resp_fr_n, pix_per_mic, 10, [2, 0.8]);         %0.22 microns per pixel for zoom = 9
if suppress_plot == 0
    figure(fig_n)
    colormap(cust_map)
    base_im = imagesc(resp_fr_n);
    hold on
    s = imagesc(diff_fr_n);
    alphamap = ((diff_fr_n - 1).*bin_pix);
    alpha(s, alphamap)
    hold off
    set(gca, 'xtick', [])
    set(gca, 'ytick', [])
    axis('square')
else
end

sig_resp_frame_n = (diff_fr_n - 1).*bin_pix;

