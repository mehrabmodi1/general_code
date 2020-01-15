%function [saved_warping_landmarks] = pick_matched_landmarks(dataset_stack, ROI_mat, tiff_times)
%This function displays an averaged image for trials 15 min apart for the user to
%click on a fixed landmark in both to map out distortions in the KC field to 
%allow ROI_mat to be warped to match it.

%Testing lines
clear all
close all
load_path = 'C:\Data\Data\Analysed_data\Suite2P_results\20190822\fly2_c739_PABAEL\1\';
reg_stack = load([load_path, 'tr_avg_stack.mat']);
dataset_stack = reg_stack.ave_stack;

ROI_mat = load([load_path, 'ROI_mat.mat']);
ROI_mat = ROI_mat.ROI_mat;
ROI_mat_flat = sum(ROI_mat, 3);
ROI_mat_flat(ROI_mat_flat > 1) = 0;     %getting rid of overlapping pixels
tiff_times = load([load_path, 'tif_time_stamps.mat']);
tiff_times = tiff_times.time_stamps;

time_window = 15;   %in minutes, the time window across which landmarks are manually matched

curr_threshm = 4.5;     %thresholds for image colormap adjustment
curr_threshm_tr1 = curr_threshm;
int_ranges = zeros(size(dataset_stack, 3), 2);


%computing list of trials time_window min apart
etime_vec = zeros( (size(tiff_times, 2) - 1), 1);     %storage vector for elapsed time between each pair of trials
for tr_n = 1:(size(tiff_times, 2) - 1)
    etime_vec(tr_n, 1) = minutes(tiff_times(tr_n + 1).tstamp - tiff_times(tr_n).tstamp);
end
etime_vec = cumsum(etime_vec);
etime_vec1 = floor(etime_vec./time_window);
[del, trn_15] = unique(etime_vec1);     %list of trials time_window minutes apart
saved_warping_landmarks.ref_tr_nums = trn_15;

trial_ni = 0;
while trial_ni < (length(trn_15) - 1)
    if exist([load_path, '\saved_warping_landmarks.mat']) == 2
        saved_warping_landmarks = load([load_path, '\saved_warping_landmarks.mat']);
        saved_warping_landmarks = saved_warping_landmarks.saved_warping_landmarks;
        threshes = saved_warping_landmarks(size(saved_warping_landmarks, 2)).threshes;
        curr_threshm_tr1 = threshes(1);
        curr_threshm_tr2 = threshes(2);
        trial_ni = size(saved_warping_landmarks, 2);
    else
    end
    
    if trial_ni >= (length(trn_15) - 1)
        break
    else
    end
    
    trial_ni = trial_ni + 1;
    trial_n1 = trn_15(trial_ni);
    trial_n2 = trn_15(trial_ni + 1);
    
    
    figure(1)
    frame_n1 = squeeze(dataset_stack(:, :, trial_n1));
    if sign(min(min(frame_n1))) == -1
        frame_n1 = frame_n1 + (-1 .* min(min(frame_n1)));
    else
    end
    [frame_obj] = plot_frame(frame_n1, curr_threshm_tr1, [1, 2, 1], []);   
    
    title('Previous trial mean and ROI. Click outside image to adjust colormap or click inside to continue.')
    done = 0;
    
    if trial_ni > 1
        curr_threshm_tr1 = curr_threshm_tr2;
    else
    end
    while done == 0
        [x1, y1] = ginput(1);
        %allowing user to ask for brighter or dimmer colormapping
        if x1 < 0 || x1 > size(frame_n1, 2) || y1 < 0 || y1 > size(frame_n1, 1)
            %pulling up options box to re-do last landmark or mark current
            %trial as z-drifted
            choice = listdlg('ListString', {'make brighter', 'make dimmer'}, 'SelectionMode', 'single');

            if choice == 1      %make brighter selected
                subplot(1, 2, 1)
                curr_threshm_tr1 = curr_threshm_tr1.*0.7;
                [frame_obj] = plot_frame(frame_n1, curr_threshm_tr1, [1, 2, 1], []);

            elseif choice == 2     %make dimmer selected
                subplot(1, 2, 1)
                curr_threshm_tr1 = curr_threshm_tr1.*1.3;
                [frame_obj] = plot_frame(frame_n1, curr_threshm_tr1, [1, 2, 1], []);

            else
            end
        else

            done = 1;
            curr_threshm = curr_threshm_tr1;
        end
    end
    
    ax = gca;
    int_ranges(trial_ni, :) = ax.CLim;


    title('Trial1 mean')

    subplot(1, 2, 2)
    frame_n2 = squeeze(dataset_stack(:, :, trial_n2));

    done = 0;
    curr_threshm_tr2 = curr_threshm_tr1;
    while done == 0
        [frame_obj] = plot_frame(frame_n2, curr_threshm_tr2, [1, 2, 2], []);
        title(['Trial ', int2str(trial_n1), ' click on image to continue or outside for more options.'  ])
        [x, y] = ginput(1);

        %checking if click was outside image
        if x < 0 || x > size(frame_n2, 2) || y < 0 || y > size(frame_n2, 1)
            %pulling up options box to re-do last landmark or mark current
            %trial as z-drifted
            choice = listdlg('ListString', {'z-drifted', 'redo landmark', 'make brighter', 'make dimmer', 'landmarks done'}, 'SelectionMode', 'single');

            if choice == 1          %z-drift selected
                z_drifted = 1;
                done = 1; 
            elseif choice == 3      %make dimmer selected
                subplot(1, 2, 2)
                curr_threshm_tr2 = curr_threshm_tr2.*0.7;
                [frame_obj] = plot_frame(frame_n2, curr_threshm_tr2, [1, 2, 2], []);
            elseif choice == 4      %make brighter selected
                subplot(1, 2, 2)
                curr_threshm_tr2 = curr_threshm_tr2.*1.3;
                [frame_obj] = plot_frame(frame_n2, curr_threshm_tr2, [1, 2, 2], []);
            elseif choice == 2      %re-do last landmark selected
                landmark_n = landmark_n - 1;
                continue
            end
        else
            
            done = 1;
        end
    end
        
    
    %asking user to specify landmark pairs 
    frames = cat(3, frame_n1, frame_n2);
    threshes = [curr_threshm_tr1; curr_threshm_tr2];
    [landmark_pairs] = capture_landmarks(frames, threshes);
    disp(['landmarks obtained for trials', int2str(trial_n1), ' and ', int2str(trial_n2), '.' ]);
    
    
    saved_warping_landmarks(trial_ni).landmark_pairs = landmark_pairs;
    saved_warping_landmarks(trial_ni).threshes = threshes;
    save([load_path, '\saved_warping_landmarks.mat'], 'saved_warping_landmarks');
    
    close figure 1

    
end
%testing lines
ROI_mat_warped = zeros(size(ROI_mat, 1), size(ROI_mat, 2));
ref_tr_nums =  saved_warping_landmarks(1).ref_tr_nums;   
for lm_setn = 1:size(saved_warping_landmarks, 2)
    curr_tr1 = ref_tr_nums(lm_setn);
    curr_im1 = dataset_stack(:, :, curr_tr1);
    curr_tr2 = ref_tr_nums(lm_setn + 1);
    curr_im2 = dataset_stack(:, :, curr_tr2);

    curr_landmark_pairs = saved_warping_landmarks(lm_setn).landmark_pairs;
    ROI_mat_warped = warp_im(ROI_mat_flat, curr_landmark_pairs);

    %plotting with overlay
    figure(1)
    plot_frame_ROI(curr_im1, threshes(1), [1, 2, 1], ROI_mat_flat);
    plot_frame_ROI(curr_im2, threshes(2), [1, 2, 2], ROI_mat_warped);
    
    %PICK UP THREAD HERE
    %figure out how to offset ROIs to match up with cells - draggable?
    keyboard
end





%This function creates a subplot with the frame and the ROI on top of it.
function [frame_obj] = plot_frame(frame, curr_thresh, subplot_n, landmark_pairs)
    figure(1)
    subplot(subplot_n(1), subplot_n(2), subplot_n(3))
    frame_obj = imagesc(frame, [0, curr_thresh.*median(reshape(frame, 1, []), 'omitnan')]);
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    colormap('gray')
    hold on
    if size(landmark_pairs, 1) > 0
        if subplot_n(3) == 1
            plot(landmark_pairs(:, 1), landmark_pairs(:, 2), 'r*')
        elseif subplot_n(3) == 2
            plot(landmark_pairs(:, 3), landmark_pairs(:, 4), 'r*')
        else
        end
    else
    end
    hold off
    plot_big_fig(1)
end


%This function creates a subplot with the frame and the ROI on top of it.
function [frame_obj] = plot_frame_ROI(frame, curr_thresh, subplot_n, ROI_mat)
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
    hold off
    plot_big_fig(1)
end


function [landmark_pairs] = capture_landmarks(frames, threshes)
    landmark_pairs = [];
    landmark_pairs_plotting = [];
    figure(1)
    all_done = 0;
    while all_done == 0
        landmark_pairs_plotting = [landmark_pairs_plotting; zeros(1, 4)];
        %getting landmark on left image
        curr_done = 0;
        curr_done2 = 0;
        while curr_done == 0
            plot_frame(frames(:, :, 1), threshes(1), [1, 2, 1], landmark_pairs);
            title('Pick landmark on left image.')
            hold on
            
            [x, y] = ginput(1);
            title(' ')
            %checking if user clicked on or outside image
            if x < 0 || x > size(frames, 2) || y < 0 || y > size(frames, 1)
                choice = listdlg('ListString', {'redo last landmark', 'Cancel', 'All done!'}, 'SelectionMode', 'single');
                
                if choice == 1
                    landmark_pairs(size(landmark_pairs, 1), :) = [];
                    landmark_pairs_plotting = landmark_pairs;
                elseif choice == 3
                    curr_done = 1;
                    curr_done2 = 1;
                    all_done = 1;
                end
                    
            else
                %case if user clicked on the left image, recording landmark
                curr_landmark_pair(1, 1) = x;
                curr_landmark_pair(1, 2) = y;
                landmark_pairs_plotting(size(landmark_pairs_plotting, 1), 1) = x;
                landmark_pairs_plotting(size(landmark_pairs_plotting, 1), 2) = y;
                curr_done = 1;
            end
            
            plot_frame(frames(:, :, 1), threshes(1), [1, 2, 1], landmark_pairs_plotting);
            
        end
        
        %getting landmark on right image
        while curr_done2 == 0
            plot_frame(frames(:, :, 2), threshes(2), [1, 2, 2], landmark_pairs);
            title('Pick matching landmark on right image.')
            hold on
            [x, y] = ginput(1);
            title(' ')
            %checking if user clicked on or outside image
            if x < 0 || x > size(frames, 2) || y < 0 || y > size(frames, 1)
                choice = listdlg('ListString', {'redo last landmark', 'Cancel', 'All done!'}, 'SelectionMode', 'single');
                if choice == 1
                    landmark_pairs(size(landmark_pairs, 1), :) = [];
                    landmark_pairs_plotting = landmark_pairs;
                elseif choice == 3
                    curr_done2 = 1;
                    all_done = 1;
                    
                                        
                end
                    
            else
                %case if user clicked on the left image, recording landmark
                curr_landmark_pair(1, 3) = x;
                curr_landmark_pair(1, 4) = y;
                landmark_pairs_plotting(size(landmark_pairs_plotting, 1), 3) = x;
                landmark_pairs_plotting(size(landmark_pairs_plotting, 1), 4) = y;
                curr_done2 = 1;
                
            end
           
            plot_frame(frames(:, :, 2), threshes(2), [1, 2, 2], landmark_pairs_plotting);
        end
        landmark_pairs = [landmark_pairs; curr_landmark_pair];
        landmark_pairs_plotting = landmark_pairs;
        
        
    end
end