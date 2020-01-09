%function [landmarks] = pick_matched_landmarks(dataset_stack, ROI_mat, tiff_times)
%This function displays an averaged image for trials 15 min apart for the user to
%click on a fixed landmark in both to map out distortions in the KC field to 
%allow ROI_mat to be warped to match it.

%Testing lines
load_path = 'C:\Data\Data\Analysed_data\Suite2P_results\20190822\fly1_c739_PABAEL\1\';
reg_stack = load([load_path, 'tr_avg_stack.mat']);
dataset_stack = reg_stack.ave_stack;

ROI_mat = load([load_path, 'ROI_mat.mat']);
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

trial_ni = 0;
while trial_ni < (length(trn_15) - 1)
    trial_ni = trial_ni + 1;
    trial_n1 = trn_15(trial_ni);
    trial_n2 = trn_15(trial_ni + 1);
    
    figure(1)
    frame_n1 = squeeze(dataset_stack(:, :, trial_n1));
    if sign(min(min(frame_n1))) == -1
        frame_n1 = frame_n1 + (-1 .* min(min(frame_n1)));
    else
    end
    [frame_obj] = plot_frame(frame_n1, curr_threshm_tr1, [1, 2, 1]);   
    
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
                curr_threshm_tr1 = curr_threshm_tr1.*0.85;
                [frame_obj] = plot_frame(frame_n1, curr_threshm_tr1, [1, 2, 1]);

            elseif choice == 2     %make dimmer selected
                subplot(1, 2, 1)
                curr_threshm_tr1 = curr_threshm_tr1.*1.15;
                [frame_obj] = plot_frame(frame_n1, curr_threshm_tr1, [1, 2, 1]);

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
        [frame_obj] = plot_frame(frame_n2, curr_threshm_tr2, [1, 2, 2]);
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
                curr_threshm_tr2 = curr_threshm_tr2.*0.85;
                [frame_obj] = plot_frame(frame_n2, curr_threshm_tr2, [1, 2, 2]);
            elseif choice == 4      %make brighter selected
                subplot(1, 2, 2)
                curr_threshm_tr2 = curr_threshm_tr2.*1.15;
                [frame_obj] = plot_frame(frame_n2, curr_threshm_tr2, [1, 2, 2]);
            elseif choice == 2      %re-do last landmark selected
                landmark_n = landmark_n - 1;
                continue
            end
        else
            
            %PICK UP THREAD HERE
            %Write function to capture landmark pairs below - check if it
            %can call plot_frame to re-plot images with landmark pairs
            %highlighted.
            
            [landmark_pairs] = capture_landmarks();


            done = 1;
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
function [frame_obj] = plot_frame(frame, curr_thresh, subplot_n)
    figure(1)
    subplot(subplot_n(1), subplot_n(2), subplot_n(3))
    frame_obj = imagesc(frame, [0, curr_thresh.*median(reshape(frame, 1, []), 'omitnan')]);
    set(gca,'xtick',[])
    set(gca,'xticklabel',[])
    set(gca,'ytick',[])
    set(gca,'yticklabel',[])
    colormap('gray')
    hold on
    plot_big_fig(1)
end

function [landmark_pairs] = capture_landmarks()
end