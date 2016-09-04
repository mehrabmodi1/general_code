%Written by Mehrab N Modi, 20140915. Calculates a difference image for
%evoked responses and identifies significantly brighter pixels in the
%difference image. Then, asks user to manually outline cell-shaped and
%sized clusters of significant pixels and defines ROIs as the filled convex
%hull of the intersection of the manually drawn ROI and the cluster of
%significantly brighter pixels.

clear all
close all


%House-keeping
%-----------------------

direc = 'C:\Data\CSHL\20141027\water_odor_stim\';                        %directory with raw data in it
save_direc = 'C:\Data\CSHL\Analysed_data\';                             %main directory to save processed data files

[no_frames, no_trials, frame_rate, zoom, d_path, f_name, tr_tg_no, ave_frame] = tiff_info(direc);
[isi, odor_t_list, odorNames, n_reps, stim_time] = load_stim_param(direc);

frame_time = 1./frame_rate;                                             %in s
stim_frame = floor(stim_time(1, 1)./frame_time);                              %frame number at which stim was delivered

zoom_cell_dia_factor = 6.9;                                             %cell dia at zoom 2.6 is typically 18 pixels (in a 256x256 pixel image)
cell_r = 0.5.*zoom_cell_dia_factor.*zoom;
cell_area = 3.14.*(cell_r.^2);

save_path = [save_direc d_path];                                        %final sub-directory to save processed data files
mkdir(save_path);

mean_max_switch = 1;                                                    %switch between taking the mean(2) or the max(2) of the post stim pixels to calculate the difference image

%loop to read in multiple trials
more_trials = 1;
if exist([save_path 'ROI_trial_n.txt']) == 2
    trial_no = load([save_path 'ROI_trial_n.txt']);
else
    trial_no = 1;
end

while more_trials == 1
    
    if trial_no == (no_trials - 1)
        more_trials = 0;
        continue
    else
    end
    
    if trial_no < (no_trials - 1)
        trial_no = trial_no + 1;
    else
    end
    
    trial_no_f = trial_no + tr_tg_no - 1;           %adding the trial tag number of the first trial of this set to get the right filename
    trial_no_f = sprintf('%03.0f', trial_no_f);
    file_path = [direc f_name trial_no_f, '.tif'];  %full file path for currently analysed trial
    
    %loop to read in all frames of current trial
    pre_frame_mat = zeros(size(ave_frame, 1), size(ave_frame, 2), (stim_frame - 1)) + nan;
    post_frame_mat = zeros(size(ave_frame, 1), size(ave_frame, 2), (no_frames - stim_frame) ) + nan;
    for frame_no = 1:no_frames
        frame = double(imread(file_path, frame_no));
        if frame_no < stim_frame
            if frame_no == 1
                pre_stim = frame;
            else
                pre_stim = pre_stim + frame;
            end
            pre_frame_mat(:, :, frame_no) = frame;
        elseif frame_no > stim_frame
            post_frame_mat(:, :, (frame_no - stim_frame) ) = frame;
        end
    end
    pre_stim = pre_stim./(stim_frame - 1);              %mean of pre-stimulus frames
    
    
    %calculating post-stimulus frame by identifying the peak value of each
    %pixel in the time series and then summing an area under the curve
    %after that time point.
    [del max_mat] = max(post_frame_mat, [], 3);
    
    
    %loop to identify area under the curve near peak response for each
    %pixel
    pre_f = round(0.1./frame_time);         %no of frames before peak to include in area calculation
    post_f = round(0.9./frame_time);        %no of frames after peak to include in area calculation
    post_pix_mat = zeros(size(frame)) + nan;
    for pix_row_n = 1:size(frame, 1)
        for pix_col_n = 1:size(frame, 2)
            pix_max_f = max_mat(pix_row_n, pix_col_n);                %frame no of max intensity for current pixel
            pre_f_curr = max([(pix_max_f - pre_f), 1]);
            post_f_curr = min([(pix_max_f + post_f), (no_frames - stim_frame)]);
            resp_vec = post_frame_mat(pix_row_n, pix_col_n, pre_f_curr:post_f_curr );        %vector of intensities on frames around max frame for current pixel
            post_pix_mat(pix_row_n, pix_col_n) = nansum(resp_vec)./length(resp_vec);
            
        end
        
    end
    
    post_stim = post_pix_mat;
    
    
    
    %identifying background pixels
    pre_stim_norm = pre_stim./max(max(pre_stim));
    bk_pix = im2bw(pre_stim_norm, 0.02);
    bk_pix = bwareaopen(bk_pix, 8);                    %getting rid of remnant, small bunches of pixels
    bk_pixi = find(bk_pix == 0);
    clear pre_stim_norm
    clear bk_pix
    
    if mean_max_switch == 1
        post_stim = post_stim(:, :, 1);                     %max value for each pix in all frames after stim frame
    elseif mean_max_switch == 2
        post_stim = post_stim(:, :, 1)./(no_frames - stim_frame);                     %mean value for each pix in all frames after stim frame
    else
    end
    
    diff_im = post_stim - pre_stim;                     %difference image, max post-stim value minus ave pre-stim baseline for each pixel
    diff_im(bk_pixi) = 0;                               %forcing background pixels to 0
    
    %loop to interactively set threshold multiplier
    try_again = 1;
    
    if exist('thresh_multiplier') == 0
        thresh_multiplier = 4;
    else
    end
    
    while try_again == 1

        thresh_mat = thresh_multiplier.*std(pre_frame_mat, 0, 3);         %threshold diff values for each pixel set as twice SD during pre-stim baseline frames
        sig_res_pixi = find(diff_im > thresh_mat);
        sig_diff_im = zeros(size(diff_im));
        sig_diff_im(sig_res_pixi) = diff_im(sig_res_pixi);                %difference image with non-significant pixels forced to 0


        %getting rid of tiny clumps of pixels
        BW = im2bw(sig_diff_im, 0.1);
        BW = bwareaopen(BW, floor(cell_area./10), 8);
        sig_diff_area = sig_diff_im.*BW; 

        %figure(1)
        scrsz = get(0,'ScreenSize');
        figure('Name','2','NumberTitle','off', 'Position', [30, 100, min(scrsz(3:4))./1.3, min(scrsz(3:4))./1.3])
        subplot(3, 1, 1)
        colormap('gray')
        imagesc(pre_stim)
        title(['averaged pre-stim frames for trial number ' int2str(trial_no)])
        
        subplot(3, 1, 2)
        colormap(gray)
        imagesc(diff_im)
        title('Difference image')
        
        subplot(3, 1, 3)
        colormap(gray)
        imagesc(sig_diff_area)
        ROI_im = sig_diff_area;
        title('Significantly brighter pixels at current threshold')
        
        %generating input dialog box
        prompt = {'Enter threshold multiplier:','If satisfied, enter 0:'};
        dlg_title = 'Manually change threshold';
        num_lines = 1;
        def = {num2str(thresh_multiplier), int2str(try_again)};
        answer = inputdlg(prompt,dlg_title,num_lines,def);
        
        try_again = str2num(answer{2, 1});
        thresh_multiplier = str2num(answer{1, 1});
        close figure 2
    end
    
    %making sure last set threh_multiplier is used
    thresh_mat = thresh_multiplier.*std(pre_frame_mat, 0, 3);         %threshold diff values for each pixel set as twice SD during pre-stim baseline frames
    sig_res_pixi = find(diff_im > thresh_mat);
    sig_diff_im = zeros(size(diff_im));
    sig_diff_im(sig_res_pixi) = diff_im(sig_res_pixi);                %difference image with non-significant pixels forced to 0
    
    
    
    %--------------------------------
    %Using manual intervention to tell cell-sized and cell-shaped objects apart

    %loop to manually draw boundaries around each cell's cluster of sig-res pixels
    more_cells = 1;
    max_val = max(max(ROI_im));
    min_val = min(min(ROI_im));
    sig_pixi = find(ROI_im > 0);
    ROI_im2 = ROI_im;
    
    while more_cells == 1
        %checking for and loading pre-existing ROI-file, or initialising a blank matrix
        if exist([save_path 'ROI_mat.txt']) == 2
            ROI_mat = load([save_path 'ROI_mat.txt']);
        else
            ROI_mat = zeros(size(ave_frame));
        end
        
        if exist([save_path 'ROI_int_mat.txt']) == 2
            ROI_int_mat = load([save_path 'ROI_int_mat.txt']);
        else
            ROI_int_mat = zeros(size(ave_frame));
        end
        
        %changing pixels in display image to indicate ROIs already chosen
        ROI_pix = find(ROI_int_mat > 0);
        ROI_im2(ROI_pix) = max_val;
        
        
        %plotting averaged pre-stim image for reference on doubtful cases
        figure(4)
        colormap('gray')
        imagesc(pre_stim)
        title('averaged pre-stim frames')
        
        %plotting image and obtaining manual ROI
        %scrsz = get(0,'ScreenSize');
        figure('Name','3','NumberTitle','off', 'Position', [10, 10, (min(scrsz(3:4)) .*1.7), min(scrsz(3:4))])
        subplot(1, 2, 1)
        imagesc(diff_im)
        title('Difference image')
        subplot(1, 2, 2)
        colormap(gray)
        imagesc(ROI_im2, [min_val, max_val])
        title('Significantly brighter pixels + ROIs')
        
        
        ROI_BW = roipoly;           %obtaining manually drawn ROI
        ROI_BW1 = ROI_BW;
        
        %finding convex hull of intersection of manually drawn ROI and significantly responsive pixels.
        ROI_pix = find(ROI_BW == 1);
        int_pixi = intersect(ROI_pix, sig_pixi);
        
        ROI_BW = zeros(size(ave_frame));
        ROI_BW(int_pixi) = 1;
        
        if sum(sum(ROI_BW)) > 0

            labels = bwlabel(ROI_BW);
            stats = regionprops(labels, 'ConvexHull');
            stats = stats.ConvexHull;
            
            x_vec = stats(:, 1);
            y_vec = stats(:, 2);
            [delete ROI_BW_f] = roifill(0, 0, ave_frame, x_vec, y_vec);         %filling in convex hull of intersection between drawn ROI and sig pixels

            ROI_int_mat = ROI_int_mat + ROI_BW_f;
            a = find(ROI_int_mat > 1);
            ROI_int_mat(a) = 0;                 %getting rid of overlapping pixels between ROIs
        else 
        end
        
        clear delete
        
        save([save_path 'ROI_mat.txt'], 'ROI_mat', '-ASCII');
        save([save_path 'ROI_int_mat.txt'], 'ROI_int_mat', '-ASCII');
        save([direc 'ROI_mat.txt'], 'ROI_int_mat', '-ASCII');
        
        
        ROI_mat = ROI_mat + ROI_BW1;            %matrix of manually drawn ROIs for saving
        
        %calling a dialog box to ask if more cells need to be drawn
        qstring = {'Mark out more cells?'};
        dlg_title = 'More cells?';
        def = {'Yes'};
        answer = questdlg(qstring,dlg_title, def);
        if strcmp(answer, 'Yes') == 1
            more_cells = 1;
        else
            more_cells = 0;
        end
        close figure 4
        close figure 3
    end
    save([save_path 'ROI_trial_n.txt'], 'trial_no', '-ASCII');
    
    %calling a dialog box to ask if more cells need to be drawn
    qstring = {'Look at more trials?'};
    dlg_title = 'More trials?';
    def = {'Yes'};
    answer = questdlg(qstring,dlg_title, def);
    if strcmp(answer, 'Yes') == 1
        more_trials = 1;
    else
        more_trials = 0;
    end
    
end


%deleting overlapping parts of ROIs manually
more_del = 1;
while more_del == 1
    figure(1)
    imagesc(ROI_int_mat)
    title('Choose area of pixels for deletion from ROIs')
    ROI_del = roipoly;
    a = find(ROI_del == 1);
    ROI_int_mat(a) = 0; 
    
    
    %calling a dialog box to ask if more cells need to be drawn
        qstring = {'Mark out more areas to delete?'};
        dlg_title = 'More deletions?';
        def = {'Yes'};
        answer = questdlg(qstring,dlg_title, def);
        if strcmp(answer, 'Yes') == 1
            more_del = 1;
        else
            more_del = 0;
        end
    
    
end
    
%saving
save([save_path 'ROI_trial_n.txt'], 'trial_no', '-ASCII');

save([save_path 'ROI_mat.txt'], 'ROI_mat', '-ASCII');
save([save_path 'ROI_int_mat.txt'], 'ROI_int_mat', '-ASCII');
save([direc 'ROI_mat.txt'], 'ROI_int_mat', '-ASCII');
 