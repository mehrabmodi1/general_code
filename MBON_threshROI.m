%direc = 'E:\Data\Analysed_data\Manual_ROIs\20200715\fly3_A1_CS_30s_trace_60s_44hrstrvd\';
direc = 'E:\Data\Analysed_data\Manual_ROIs\20200716\fly3_G2_PABAEL_handover_noLED_nodendrites\';
ave_stack = load([direc, 'tr_avg_stack.mat']);
ave_stack = ave_stack.ave_stack;
curr_thresh = 4.5;
for trial_n = 1:size(ave_stack, 3)
    curr_im = squeeze(ave_stack(:, :, trial_n));
    curr_im = curr_im./max(max(curr_im));
    
    figure(1)
    %plotting image alone
    subplot(1, 2, 1)
    plot_frame(curr_im, curr_thresh, [1, 2, 1], zeros(size(curr_im, 1), size(curr_im, 2)));
    
    
    %filtering image
    se = strel('disk', round(size(curr_im, 1)./50));
    curr_im_f1 = imtophat(curr_im, se);
    curr_im_f2 = imgaussfilt(curr_im_f1, 8);
    
%     figure(2)
%     imagesc(curr_im_f2)
    %thresholding and identifying ROI
    T = graythresh(curr_im_f2);
    ROI = imbinarize(curr_im_f2, T);
    %plotting image with ROI
    plot_frame(curr_im, curr_thresh, [1, 2, 2], ROI);
    
    if rem(trial_n, 5) == 0

        keyboard
    else
    end
%     %plotting results of currently detected ROI
%     plot_frame(curr_im, curr_thresh, [], ROI);
    
end

%This function creates a subplot with the frame and the ROI on top of it.
function [frame_obj, ROI_obj] = plot_frame(frame, curr_thresh, subplot_n, ROI_mat)
    
    if isempty(subplot_n) ~= 1
        fig_n = 1;
        figure(fig_n)
        subplot(subplot_n(1), subplot_n(2), subplot_n(3))
    else
        fig_n = 2;
    end
    
    figure(fig_n)
    if median(reshape(frame, 1, []), 'omitnan') ~= 0
        frame_obj = imagesc(frame, [0, curr_thresh.*median(reshape(frame, 1, []), 'omitnan')]);
    elseif median(reshape(frame, 1, []), 'omitnan') == 0
        frame_obj = imagesc(frame, [0, 1]);
    end
    
    
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    colormap('gray')
    hold on
    ROI_mat_sc = ROI_mat.*max(max(frame));
    ROI_obj = imagesc(ROI_mat_sc);
    hold off
    ROI_obj.AlphaData = ROI_mat.*0.5;
    plot_big_fig(fig_n)
end