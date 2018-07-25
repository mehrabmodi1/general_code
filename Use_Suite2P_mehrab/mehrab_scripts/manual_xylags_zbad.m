function [lag_mat, bad_trs] = manual_xylags_zbad(dataset_stack)
%This function displays an averaged image of each frame for the user to
%click on a fixed landmark to determine x-y lags for slow drift as well as
%to allow the user to report bad, z-drift trials by clicking outside the
%image.

%displaying averaged frame
lag_mat = zeros(size(dataset_stack, 3), 2);
bad_trs = zeros(size(dataset_stack, 3), 1);

for trial_n = 1:size(dataset_stack, 3)

    figure(1)
    subplot(1, 2, 1)
    frame1 = squeeze(dataset_stack(:, :, 1));
    imagesc(frame1, [0, 0.2.*max(max(frame1))]);
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    colormap('gray')
    plot_big_fig(1);
    
    if trial_n == 1
        title('Trial1 mean, click on a landmark.')
        [x1, y1] = ginputc(1, 'Color', [1, 0, 0]);
        landmark_ROI = zeros(size(frame1, 1), size(frame1, 2));
        landmark_ROI = draw_circle(y1, x1, 3, landmark_ROI, 1);
        landmark_ROI_rgb = landmark_ROI;
        landmark_ROI_rgb(:, :, 2:3) = zeros(size(frame1, 1), size(frame1, 2), 2);
    else
        title('Trial1 mean with landmark highlighted.')
        hold on
        h = imagesc(landmark_ROI_rgb);
        h.AlphaData = landmark_ROI;
        hold off
        
        subplot(1, 2, 2)
        curr_frame = squeeze(dataset_stack(:, :, trial_n));
        imagesc(curr_frame, [0, 0.2.*max(max(curr_frame))])
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
                choice = questdlg('What would you like to do?', 'landmark-marking', 'z-drifted', 'redo landmark', 'z-drifted');
                if strcmp(choice, 'z-drifted') == 1
                    z_drifted = 1;
                    done = 1;
                else
                end
            else
                z_drifted = 0;
                done = 1;
            end
        end
        if z_drifted == 0
            lag_mat(trial_n, 1) = x - x1;
            lag_mat(trial_n, 2) = y - y1;
            lag_mat(trial_n, 3) = 0;
            z_drifted = 0;
        elseif z_drifted == 1
            lag_mat(trial_n, 3) = 1;
        end
        
    end
   
end
close figure 1
