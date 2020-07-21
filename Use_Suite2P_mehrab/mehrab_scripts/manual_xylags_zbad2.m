function [lag_mat, bad_trs, done_marking, bk_ROI, ROI_mat_adj] = manual_xylags_zbad2(dataset_stack, ROI_mat, ROI_mat_warped, warping_landmarks, ROI_dilate)
%This function displays an averaged image of each frame for the user to
%click on a fixed landmark to determine x-y lags for slow drift as well as
%to allow the user to report bad, z-drift trials by clicking outside the
%image.
lag_mat = zeros(size(dataset_stack, 3), 3);
bad_trs = zeros(size(dataset_stack, 3), 1);
reg_stack = zeros(size(dataset_stack, 1), size(dataset_stack, 2), size(dataset_stack, 3));
curr_threshm = 4.5;
curr_threshm_tr1 = curr_threshm;
int_ranges = zeros(size(dataset_stack, 3), 2);

if isempty(warping_landmarks) == 1 && size(ROI_mat, 3) ~= 1
    warped_ROIs = 0;
    ROI_mat_adj = [];   %only applicable to single, large, MBON-like ROIs
elseif isempty(warping_landmarks) == 0
    warped_ROIs = 1;
    ROI_mat_adj = [];   %only applicable for single, large, MBON-like ROIs 
elseif size(ROI_mat, 3) == 1        
    warped_ROIs = 2;        %case when ROIs are of a single MBON
    disp('WARNING: current dataset auto-detected as an MBON dataset.')
    ROI_mat_adj = repmat(ROI_mat, 1, 1, size(dataset_stack, 3));      %creating an ROI matrix with one ROI for each trial and initializing to trial 1 ROI
    
else
end


trial_n = 0;
global im_posx1
global im_posy1
last_lagx = 0;
last_lagy = 0;

while trial_n < size(dataset_stack, 3)
    trial_n = trial_n + 1;
    figure(1)
    frame1 = squeeze(dataset_stack(:, :, 1));
    if sign(min(min(frame1))) == -1
        frame1 = frame1 + (-1 .* min(min(frame1)));
    else
    end
    
    %swapping in appropriate, warped ROI_mat for the original ROI_mat
    if warped_ROIs == 1
        ref_trs = warping_landmarks(1).ref_tr_nums;     %list of waypoint trials 15 min apart for which landmarks were picked
        wpnts_passed = sum(trial_n > ref_trs);
        if trial_n == 1
            wpnts_passed = 1;
        elseif wpnts_passed == length(ref_trs)      %case when a few more trials exist past the last waypoint.
            wpnts_passed = wpnts_passed - 1;
        end
        
        try
            ROI_mat = ROI_mat_warped(:, :, wpnts_passed);
        catch
            keyboard
        end
    
    elseif warped_ROIs == 0
    elseif warped_ROIs == 2
        ROI_mat = ROI_mat_adj(:, :, trial_n);
    else
    end
    
    
        
    disp('Slow, manual motion correction beginning...')
    if trial_n == 1
        [frame_obj, ROI_obj] = plot_frame(frame1, curr_threshm, [1, 2, 1], ROI_mat);      
        %recording original position of frame 1 image displayed.
        im_posx0 = ROI_obj.XData;
        im_posx0 = im_posx0(1);
        im_posy0 = ROI_obj.YData;
        im_posy0 = im_posy0(1);
        
         
        title('Draw background ROI.')
        bk_ROI = roipoly;
                
        title('Trial1 mean and ROI. Click to continue.')
        done = 0;
        
        while done == 0
            [x1, y1] = ginput(1);
            %allowing user to ask for brighter or dimmer colormapping
            if x1 < 0 || x1 > size(frame1, 2) || y1 < 0 || y1 > size(frame1, 1)
                %pulling up options box to re-do last landmark or mark current
                %trial as z-drifted
                choice = listdlg('ListString', {'make brighter', 'make dimmer'}, 'SelectionMode', 'single');
                ROI_mat_tr1 = ROI_mat;
                if choice == 1      %make brighter selected
                    subplot(1, 2, 1)
                    curr_threshm_tr1 = curr_threshm_tr1.*0.85;
                    [frame_obj, ROI_obj] = plot_frame(frame1, curr_threshm_tr1, [1, 2, 1], ROI_mat_tr1);
                    
                elseif choice == 2     %make dimmer selected
                    subplot(1, 2, 1)
                    curr_threshm_tr1 = curr_threshm_tr1.*1.15;
                    [frame_obj, ROI_obj] = plot_frame(frame1, curr_threshm_tr1, [1, 2, 1], ROI_mat_tr1);
                    
                else
                end
            else
                
                x1 = im_posx0;
                y1 = im_posy0;
                done = 1;
                curr_threshm = curr_threshm_tr1;
            end
        end
        
        reg_stack(:, :, 1) = frame1;        %first frame of registered stack is the reference frame from averaging trial 1
        ax = gca;
        int_ranges(trial_n, :) = ax.CLim;
        
    else   %Now dealing with trial_n > 1
        title('Trial1 mean with ROI.')
                
        subplot(1, 2, 2)
        curr_frame_orig = squeeze(dataset_stack(:, :, trial_n));
        curr_frame = translate_stack (curr_frame_orig, [last_lagy; last_lagx], nan);
        
        done = 0;
        while done == 0
            try
                [frame_obj, ROI_obj] = plot_frame(curr_frame, curr_threshm, [1, 2, 2], ROI_mat);
            catch
                keyboard
            end
            title(['Trial ', int2str(trial_n), ' mean, drag to match ROI, or click to bring up cursor.'  ])
            draggable(ROI_obj, 'none', [-inf inf -inf inf], 'endfcn', @end_drag_func); 
            uiwait(gcf)
            title(['Trial ', int2str(trial_n), ' click on image to continue or outside for more options.'  ])
            [x, y] = ginput(1);
            
            %checking if click was outside image
            if x < 0 || x > size(frame1, 2) || y < 0 || y > size(frame1, 1)
                %pulling up options box to re-do last landmark or mark current
                %trial as z-drifted
               
                if warped_ROIs == 2 %case where an MBON dataset was detected, allowing manual ROI adjustment due to tissue warping
                    choice = listdlg('ListString', {'z-drifted', 'redo landmark', 'make brighter', 'make dimmer', 'adjust ROI'}, 'SelectionMode', 'single');
                else
                    choice = listdlg('ListString', {'z-drifted', 'redo landmark', 'make brighter', 'make dimmer'}, 'SelectionMode', 'single');
                end
                    
                if choice == 1          %z-drift selected
                    z_drifted = 1;
                    done = 1; 
                elseif choice == 3      %make brighter selected
                    subplot(1, 2, 2)
                    curr_threshm = curr_threshm.*0.85;
                    [frame_obj, ROI_obj] = plot_frame(curr_frame, curr_threshm, [1, 2, 2], ROI_mat);
                elseif choice == 4      %make dimmer selected
                    subplot(1, 2, 2)
                    curr_threshm = curr_threshm.*1.15;
                    [frame_obj, ROI_obj] = plot_frame(curr_frame, curr_threshm, [1, 2, 2], ROI_mat);
                elseif choice == 2      %re-do last landmark selected
                    trial_n = trial_n - 1;
                    continue
                elseif choice == 5      %adjust ROI selected - used to manually adjust single/MBON ROIs
                   
                   
                    figure(2)
                    plot_frame(curr_frame, curr_threshm, [], zeros(size(ROI_mat, 1), size(ROI_mat, 2)) );
                    curr_ROI_adj = roipoly();
                    if ROI_dilate > 0
                        str = strel('disk', ROI_dilate, 0); 
                        curr_ROI_adj = imdilate(curr_ROI_adj, str);
                    else
                    end
                    ROI_mat = curr_ROI_adj;
                    ROI_mat_adj(:, :, trial_n:end) = repmat(curr_ROI_adj, 1, 1, (size(dataset_stack, 3) - (trial_n - 1) ));
                    close figure 2
                else
                end
                
                   
                
            else
                
                z_drifted = 0;
                x = im_posx1;
                y = im_posy1;
                x_lag = im_posx0 - im_posx1 + last_lagx;        %adding last tr's lag because this is automatically applied without any dragging
                y_lag = im_posy0 - im_posy1 + last_lagy;        %adding last tr's lag because this is automatically applied without any dragging
                
                %keeping track of last lags to plot next image
                last_lagx = x_lag;
                last_lagy = y_lag;
                
                curr_frame_reg = translate_stack (curr_frame_orig, [y_lag; x_lag], nan);
                
                done = 1;
       
            end
        end
        if z_drifted == 0
            lag_mat(trial_n, 1) = y_lag;
            lag_mat(trial_n, 2) = x_lag;
            lag_mat(trial_n, 3) = 0;
            z_drifted = 0;
            
            %saving a lag-corrected version of the current frame
            reg_stack(:, :, (trial_n - sum(bad_trs))) = curr_frame_reg;
            ax = gca;
            int_ranges(trial_n, :) = ax.CLim;
           
        elseif z_drifted == 1
            lag_mat(trial_n, 3) = 1;
            bad_trs(trial_n, 1) = 1;
            %int_ranges(trial_n, :) = [];
        end
        
    end
    
    
end
bad_tr_list = find(bad_trs == 1);
reg_stack(:, :, bad_tr_list) = [];

close figure 1

%playing back trial frames for manual review
playStack_specint(reg_stack, 30, int_ranges)

choice = questdlg('Alignment OK?', 'Reviewing alignment', 'Yes, stack is OK', 'redo landmarks', 'Yes, stack is OK');
if strcmp(choice, 'Yes, stack is OK') == 1
    done_marking = 1;
elseif strcmp(choice, 'redo landmarks') == 1
    done_marking = 0;
else
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
    
    
%This function is called by draggable every time the user drags and then
%drops the ROI after translating it
function end_drag_func(ROI_obj)
    global im_posx1
    global im_posy1
    im_posx1 = ROI_obj.XData;
    im_posx1 = im_posx1(1);
    im_posy1 = ROI_obj.YData;
    im_posy1 = im_posy1(1); 
    uiresume(gcf)
