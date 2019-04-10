clear all
close all

list_direc = ['D:\Data\CSHL\dataset_list_fluorescein_diffusion_20150803.txt'];
fid = fopen(list_direc);
direc_counter = 0;

%loop to go through all experiment datasets listed in list file
while 1
    
    direc = fgetl(fid);
       
    if ischar(direc) ~= 1
        break
    else
    end
    
    direc = [direc '\'];
    
    %replacing C: with D:
    a = findstr(direc, 'C:');
    direc(a) = 'D';
    
    cd(direc)
    dir_contents = dir('*.tif');
    
    direc_counter = direc_counter + 1;
    
    n_trials = size(dir_contents, 1);
    
    f_name = dir_contents(1).name;
    info = imfinfo(f_name);
    info = info(1).ImageDescription;
    n_framesi = findstr(info, 'state.acq.numberOfFrames');
    n_frames = info((n_framesi+ 25):(n_framesi + 27));
    n_frames = str2num(n_frames);
       
    %reading in data from each trial
    data_mat = zeros(5, n_trials);
    raw_resp_traces = zeros(n_frames, n_trials) + nan;
    
    try 
        ROI_mat = load('ROI_mat.txt');
        ROIs_saved = 1;
    catch
        ROIs_saved = 0;
    end
    
    for trial_n = 1:n_trials
       f_name =  dir_contents(trial_n).name;
       
       %asking user to draw ROIs around background patches
       %averaged image obtaining
       if trial_n == 1 && ROIs_saved == 0
           for frame_n = 1:n_frames;
               if frame_n == 1
                    ave_frame = double(imread(f_name, frame_n));
               elseif frame_n > 1
                   ave_frame = ave_frame + double(imread(f_name, frame_n));
               else
               end
               
           end
           ave_frame = ave_frame./n_frames;
           
           figure(1)
           imagesc(ave_frame)
           colormap('gray');
           ROI_mat = zeros(size(ave_frame, 1), size(ave_frame, 2));
           %drawing background and cell bod + calyx ROIs
           for bk_roi_n = 1:4
               if bk_roi_n < 4
                   disp('Draw Background ROI')
               else
                   disp('Draw cell body + calyx ROI')
               end
               BW = roipoly();
               ROI_mat(:, :) = ROI_mat(:, :) + (BW.*bk_roi_n);
           end
           
           save('ROI_mat.txt', 'ROI_mat', '-ASCII');           
       else
       end
       
       
       
       %reading in raw fluorescence data
       data_vec = zeros(5, n_frames) + nan;
       for frame_n = 1:n_frames
           frame = imread(f_name, frame_n);
           
           for curr_roi_n = 1:4
               curr_roi_i = find(ROI_mat == curr_roi_n);
               data_vec(curr_roi_n, frame_n) = mean(mean(frame(curr_roi_i)));
               
               if curr_roi_n == 4
                   raw_resp_traces(frame_n, trial_n) = data_vec(curr_roi_n, frame_n);
               else
               end
               
           end

       end
       data_vec = mean(data_vec, 2);
       
       
       info_t = imfinfo(f_name);
       info_t = info_t(1).ImageDescription;
       t_i = findstr(info_t, 'triggerTimeString');
       tr_time_curr = info_t((t_i + 29):(t_i + 40));
       if trial_n == 1
           tr_time1i = tr_time_curr;
           tr_time1 = datevec(tr_time_curr);
           data_vec(5) = 0;
       elseif trial_n > 1
           t_diff = etime(datevec(tr_time_curr), tr_time1);
           data_vec(5) = t_diff;
           
       else
       end
       
       data_mat(:, trial_n) = data_vec;
       
    end
    
    
    %raw plots
    figure(2)
    plot(data_mat(5, :), data_mat(1:4, :), 'LineWidth', 3)
    hold on
    
    
    %Normalised plots, bk only
    figure(3)
    max_vals = max(data_mat, [], 2);
    max_vals(5) = 1;                                             %don't want to normalise time-point info
    max_vals_r = repmat(max_vals, 1, size(data_mat, 2) );
    norm_data_mat = data_mat./max_vals_r;                        %normalising raw fluorescence data
    plot(norm_data_mat(5, :), norm_data_mat(1:3, :))
    hold on
    
    %normalised plots, 1st bk ROI only
    figure(4)
    
    %getting rid of a bad trial
    if direc_counter == 6
        norm_data_mat(:, end) = [];
        
    else
    end
    
    plot(norm_data_mat(5, :), norm_data_mat([1,3], :), 'LineWidth', 3)
    xlabel('Time (s)')
    ylabel('Norm. fluorescence intensity near KCs')
    hold on
    
    
    %fitting exponentials to curves
    curve_data = nanmean(norm_data_mat([1,3], :), 1);
    curve_data_f = 1 - curve_data;
    fitobj = fit(norm_data_mat(5, :)', curve_data_f', 'exp1');
    
    a = fitobj.a;       %coefficients in the eqn y = a*exp(b*x)
    b = fitobj.b;
    evaluated_points = a.*exp(b.*norm_data_mat(5, :));
    evaluated_points = 1 - evaluated_points;
    
    figure(5)
    plot(norm_data_mat(5, :), curve_data, '.', 'MarkerSize', 15)
    hold on
    plot(norm_data_mat(5, :), evaluated_points, 'LineWidth', 3)
    xlabel('time (s)')
    ylabel('normalized fluorescence')
end
fclose(fid);