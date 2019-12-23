%function [landmarks] = pick_matched_landmarks(dataset_stack, ROI_mat, tiff_times)
%This function displays an averaged image for trials 15 min apart for the user to
%click on a fixed landmark in both to map out distortions in the KC field to 
%allow ROI_mat to be warped to match it.

%Testing lines
load_path = 'C:\Data\Data\Analysed_data\Suite2P_results\20190822\fly1_c739_PABAEL\1\';
dataset_stack = load([load_path, 'tr_avg_stack.mat']);

ROI_mat = load([load_path, 'ROI_mat.mat']);
tiff_times = load([load_path, 'tif_time_stamps.mat']);
keyboard

curr_threshm = 4.5;
curr_threshm_tr1 = curr_threshm;
int_ranges = zeros(size(dataset_stack, 3), 2);
global im_posx1
global im_posy1

%PICK UP THREAD HERE
%computing list of trials 15 min apart


while trial_n < size(dataset_stack, 3)
    trial_n = trial_n + 1;
    figure(1)
    frame1 = squeeze(dataset_stack(:, :, 1));
    if sign(min(min(frame1))) == -1
        frame1 = frame1 + (-1 .* min(min(frame1)));
    else
    end
    
    [frame_obj, ROI_obj] = plot_frame(frame1, curr_threshm, [1, 2, 1], ROI_mat);   
    %recording original position of frame 1 image displayed.
    im_posx0 = ROI_obj.XData;
    im_posx0 = im_posx0(1);
    im_posy0 = ROI_obj.YData;
    im_posy0 = im_posy0(1);
        
    disp('Slow, manual motion correction beginning...')
    if trial_n == 1
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
                
                if choice == 1      %make brighter selected
                    subplot(1, 2, 1)
                    curr_threshm_tr1 = curr_threshm_tr1.*0.85;
                    [frame_obj, ROI_obj] = plot_frame(frame1, curr_threshm_tr1, [1, 2, 1], ROI_mat);
                    
                elseif choice == 2     %make dimmer selected
                    subplot(1, 2, 1)
                    curr_threshm_tr1 = curr_threshm_tr1.*1.15;
                    [frame_obj, ROI_obj] = plot_frame(frame1, curr_threshm_tr1, [1, 2, 1], ROI_mat);
                    
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
            [frame_obj, ROI_obj] = plot_frame(curr_frame, curr_threshm, [1, 2, 2], ROI_mat);
            title(['Trial ', int2str(trial_n), ' mean, drag to match ROI, or click to bring up cursor.'  ])
            draggable(ROI_obj, 'none', [-inf inf -inf inf], 'endfcn', @end_drag_func); 
            uiwait(gcf)
            title(['Trial ', int2str(trial_n), ' click on image to continue or outside for more options.'  ])
            [x, y] = ginput(1);
            
            %checking if click was outside image
            if x < 0 || x > size(frame1, 2) || y < 0 || y > size(frame1, 1)
                %pulling up options box to re-do last landmark or mark current
                %trial as z-drifted
                choice = listdlg('ListString', {'z-drifted', 'redo landmark', 'make brighter', 'make dimmer'}, 'SelectionMode', 'single');
                
                if choice == 1          %z-drift selected
                    z_drifted = 1;
                    done = 1; 
                elseif choice == 3      %make dimmer selected
                    subplot(1, 2, 2)
                    curr_threshm = curr_threshm.*0.85;
                    [frame_obj, ROI_obj] = plot_frame(curr_frame, curr_threshm, [1, 2, 2], ROI_mat);
                elseif choice == 4      %make brighter selected
                    subplot(1, 2, 2)
                    curr_threshm = curr_threshm.*1.15;
                    [frame_obj, ROI_obj] = plot_frame(curr_frame, curr_threshm, [1, 2, 2], ROI_mat);
                elseif choice == 2      %re-do last landmark selected
                    trial_n = trial_n - 1;
                    continue
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
    figure(1)
    subplot(subplot_n(1), subplot_n(2), subplot_n(3))
    frame_obj = imagesc(frame, [0, curr_thresh.*median(reshape(frame, 1, []), 'omitnan')]);
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
    plot_big_fig(1)
end
    
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
end
