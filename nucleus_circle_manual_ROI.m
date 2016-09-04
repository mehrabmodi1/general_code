%Written by Mehrab N Modi, 20141110. Calculates an averaged image for the
%first trial, asks user to mark centers of cell nuclei and draws circluar
%ROIs at each of thse points.

clear all
close all

list_path = 'C:\Data\CSHL\dataset_list_20141121.txt';
saving = 1;

fid=fopen(list_path);
dir_counter = 0;
while 1
    dir_counter = dir_counter + 1;
    direc = fgetl(fid);
    if ~ischar(direc),   break,   end
    direc = [direc '\'];
    
    
    disp(['dataset number ' int2str(dir_counter)])


    %House-keeping
    %-----------------------
    save_direc = 'C:\Data\CSHL\Analysed_data\';                             %main directory to save processed data files

    [no_frames, no_trials, frame_rate, zoom, d_path, f_name, tr_tg_num, ave_frame] = tiff_info(direc);
    [isi, odor_t_list, odorNames, n_reps, stim_time] = load_stim_param(direc);

    frame_time = 1./frame_rate;                                             %in s
    stim_frame = floor(stim_time(1, 1)./frame_time);                              %frame number at which stim was delivered

    zoom_cell_dia_factor = 6.4;                                             %cell dia at zoom 2.6 is typically 18 pixels (in a 256x256 pixel image)
    cell_r = 0.5.*zoom_cell_dia_factor.*zoom;
    cell_area = 3.14.*(cell_r.^2);

    save_path = [save_direc d_path];                                        %final sub-directory to save processed data files
    mkdir(save_path);

    %obtaining high-quality averaged image from first trial
    trial_n = 1;
    trial_no_f = trial_n + tr_tg_num - 1;           %adding the trial tag number of the first trial of this set to get the right filename
    trial_no_f = sprintf('%03.0f', trial_no_f);
    file_path = [direc f_name trial_no_f, '.tif'];  %full file path for currently analysed trial
    counter = 0;
    thresh = 0.7;
    c_vec = zeros(1, 30);
    for frame_n = 5:35
        frame = double(imread([file_path], frame_n));

        if frame_n == 5
            ave_frame = frame;
            frame1 = frame;
            continue
        else
        end

        %aligning current frame to frames read in so far
        c = normxcorr2(ave_frame, frame1);

        if max(max(c)) < thresh
            continue
        else
        end

        frame = xcorr2_aligner(frame, c);

        c_vec(frame_n) = max(max(c));

        ave_frame = ave_frame + frame;
        counter = counter + 1;
        
    end
    ave_frame = ave_frame./counter;

    
    %manually obtaining ROIs
    figure(1)
    colormap('gray')
    pix_min = min(min(ave_frame));
    pix_max = max(max(ave_frame));
    imagesc(ave_frame, [pix_min, (pix_max.*1.4)])
    set(gcf, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);

    figure(2)
    plot(c_vec)
    
    
    %checking for and loading pre-existing ROI-file, or initialising a blank matrix
    if exist([save_path 'ROI_mat_man.txt']) == 2
        ROI_mat = load([save_path 'ROI_mat_man.txt']);
    else
        ROI_mat = zeros(size(ave_frame, 1), size(ave_frame, 2));
    end
    
    
    cell_counter = 0;
    n_skip_to = 1;
    more_cells = 1;
    ave_frame2 = ave_frame;
    while more_cells == 1
        cell_counter = cell_counter + 1;
               
        %changing pixels in display image to indicate ROIs already chosen
        ROI_pix = find(ROI_mat > 0);
        ave_frame2(ROI_pix) = pix_max.*1.4; 
        
        figure(1)
        colormap('gray')
        imagesc(ave_frame2, [pix_min, (pix_max.*1.4)])
        set(gcf, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);
        
        %obtaining mouse pointer position on click
        [y, x] = ginput(1);
        
        curr_cell_ROI = draw_circle(zeros(size(ave_frame, 1), size(ave_frame, 2)), x, y, cell_r./2.3);
        
        %checking if adding a circle of 1 cell dia overlaps with other ROIs
        %and drawing a smaller circle if this is true.
        ROI_mat2 = ROI_mat + curr_cell_ROI;
        
        
        %getting rid of overlapping portion of ROIs, plus some more at the
        %edges
        if max(max(ROI_mat2)) > 1
            temp = find(ROI_mat2 > 1);
            overlap_mat = zeros(size(ROI_mat));
            overlap_mat(temp) = 1;
            overlap_mat = imdilate(overlap_mat, strel('disk', 3) );
            
            %stretching patches of overlap along their major axis to get
            %rid of connecting pixels
            overlap_mat1 = stretch_maj_axis(overlap_mat, 8);
            
            temp = find(overlap_mat1 > 0);
            
            %forcing overlap patch pixels to 0 in ROI mat
            ROI_mat2(temp) = 0;
        else
        end
        
        ROI_mat = ROI_mat2;
        
                       
        if cell_counter == n_skip_to
            button = questdlg('More cells?','Select more cells?','~25 more','1 more', 'done', 'done');
        
        else
        end
        
        %parsing GUI input
        if strcmp(button, '~25 more') == 1
            n_skip_to = cell_counter + 25;
            button = 'a';
        elseif strcmp(button, '1 more') == 1
            n_skip_to = cell_counter + 1;
            button = 'a';
        elseif strcmp(button, 'done')
            more_cells = 0;
        else
        end
        
        if saving == 1
            %saving ROI mat to file
            save([save_path 'ROI_mat_man.txt'], 'ROI_mat', '-ASCII');
        else
        end

    end
    
    %quality control - allows user to manually remove pixels from ROI
    %matrix
    more_edits = 1;
    while more_edits == 1

        ave_frame2 = ave_frame;
        temp = find(ROI_mat > 0);
        ave_frame2(temp) = pix_max;
        

        figure(2)
        colormap('gray')
        imagesc(ave_frame2, [pix_min, (pix_max.*1.4)])
        set(gcf, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);
        
        BW = roipoly;
        temp = find(BW == 1);
        
        ROI_mat(temp) = 0;
        
        button = questdlg('More ROI edits?')
        
        if strcmp(button, 'No') == 1 | strcmp(button, 'Cancel') == 1
            more_edits = 0;
        else
        end
        
    end
        
    
    
end
fclose(fid);