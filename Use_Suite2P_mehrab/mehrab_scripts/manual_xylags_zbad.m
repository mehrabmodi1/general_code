function [lag_mat, bad_trs, done_marking, bk_ROI] = manual_xylags_zbad(dataset_stack)
%This function displays an averaged image of each frame for the user to
%click on a fixed landmark to determine x-y lags for slow drift as well as
%to allow the user to report bad, z-drift trials by clicking outside the
%image.

lag_mat = zeros(size(dataset_stack, 3), 3);
bad_trs = zeros(size(dataset_stack, 3), 1);
reg_stack = zeros(size(dataset_stack, 1), size(dataset_stack, 2), size(dataset_stack, 3));
curr_threshm = 2;
curr_threshm_tr1 = curr_threshm;
int_ranges = zeros(size(dataset_stack, 3), 2);
trial_n = 0;
while trial_n < size(dataset_stack, 3)
    trial_n = trial_n + 1;
    figure(1)
    frame1 = squeeze(dataset_stack(:, :, 1));
    if sign(min(min(frame1))) == -1
        frame1 = frame1 + (-1 .* min(min(frame1)));
    else
    end
    
    plot_frame(frame1, curr_threshm, [1, 2, 1])   
    disp('Slow, manual motion correction beginning...')
    if trial_n == 1
        title('Draw background ROI.')
        bk_ROI = roipoly;
                
        title('Trial1 mean, click on a landmark.')
        done = 0;
       
        while done == 0
            [x1, y1] = ginputc(1, 'Color', [1, 0, 0]);
            %allowing user to ask for brighter or dimmer colormapping
            if x1 < 0 || x1 > size(frame1, 2) || y1 < 0 || y1 > size(frame1, 1)
                %pulling up options box to re-do last landmark or mark current
                %trial as z-drifted
                choice = listdlg('ListString', {'make brighter', 'make dimmer'}, 'SelectionMode', 'single');

                if choice == 1      %make dimmer selected
                    subplot(1, 2, 1)
                    curr_threshm_tr1 = curr_threshm_tr1.*0.85;
                    imagesc(frame1, [0, curr_threshm_tr1.*median(reshape(frame1, 1, []))])
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                    set(gca,'ytick',[])
                    set(gca,'yticklabel',[])
                elseif choice == 2     %make brighter selected
                    subplot(1, 2, 1)
                    curr_threshm_tr1 = curr_threshm_tr1.*1.15;
                    imagesc(frame1, [0, curr_threshm_tr1.*median(reshape(frame1, 1, []))])
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                    set(gca,'ytick',[])
                    set(gca,'yticklabel',[])
                else
                end
            else
                done = 1;
            end
        end
        landmark_ROI = zeros(size(frame1, 1), size(frame1, 2));
        landmark_ROI = draw_circle(y1, x1, 3, landmark_ROI, 1);
        landmark_ROI_rgb = landmark_ROI;
        landmark_ROI_rgb(:, :, 2:3) = zeros(size(frame1, 1), size(frame1, 2), 2);
        reg_stack(:, :, 1) = frame1;
        ax = gca;
        int_ranges(trial_n, :) = ax.CLim;
        
    else
        title('Trial1 mean with landmark highlighted.')
        hold on
        h = imagesc(landmark_ROI_rgb);
        h.AlphaData = landmark_ROI;
        hold off
        
        subplot(1, 2, 2)
        curr_frame = squeeze(dataset_stack(:, :, trial_n));
        %curr_threshm = 0.05;
        imagesc(curr_frame, [0, curr_threshm.*median(reshape(curr_frame, 1, []))])
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        colormap('gray')
        title(['Trial ', int2str(trial_n), ' mean, click on landmark, or outside image if z-drifted.'  ])
                
        %checking if click was outside image
        done = 0;
        while done == 0
            [x, y] = ginputc(1, 'Color', [1, 0, 0]);
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
                    imagesc(curr_frame, [0, curr_threshm.*median(reshape(curr_frame, 1, []))])
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                    set(gca,'ytick',[])
                    set(gca,'yticklabel',[])
                elseif choice == 4      %make brighter selected
                    subplot(1, 2, 2)
                    curr_threshm = curr_threshm.*1.15;
                    imagesc(curr_frame, [0, curr_threshm.*median(reshape(curr_frame, 1, []))])
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                    set(gca,'ytick',[])
                    set(gca,'yticklabel',[])
                elseif choice == 2      %re-do last landmark selected
                    trial_n = trial_n - 1;
                    continue
                end
            else
                z_drifted = 0;
                done = 1;
            end
        end
        if z_drifted == 0
            lag_mat(trial_n, 1) = y1 - y;
            lag_mat(trial_n, 2) = x1 - x;
            lag_mat(trial_n, 3) = 0;
            z_drifted = 0;
            %saving a lag-corrected version of the current frame
            reg_stack(:, :, (trial_n - sum(bad_trs))) = translate_stack(curr_frame, [lag_mat(trial_n, 1); lag_mat(trial_n, 2)], nan);
            ax = gca;
            int_ranges(trial_n, :) = ax.CLim;
           
        elseif z_drifted == 1
            lag_mat(trial_n, 3) = 1;
            bad_trs(trial_n, 1) = 1;
            int_ranges(trial_n, :) = [];
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

function plot_frame(frame, curr_thresh, subplot_n)
    figure(1)
    subplot(subplot_n(1), subplot_n(2), subplot_n(3))
    imagesc(frame, [0, curr_thresh.*median(reshape(frame, 1, []))])
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    colormap('gray')
    plot_big_fig(1)