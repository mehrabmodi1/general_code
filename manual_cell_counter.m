clear all
close all

list_path = 'C:\Data\CSHL\dataset_list_20141031.txt';
save_direc = 'C:\Data\CSHL\Analysed_data\';                             %main directory to save processed data files

fid=fopen(list_path);
dir_counter = 0;
while 1
    dir_counter = dir_counter + 1;
    direc = fgetl(fid);
    if ~ischar(direc),   break,   end
    direc = [direc '\'];
    
    
    disp(['dataset number ' int2str(dir_counter)])
    
    
    %house-keeping
    [n_frames, n_trials, frame_rate, zoom, d_path, f_name, tr_tg_num, ave_frame] = tiff_info(direc);
    save_path = [save_direc d_path];

    trial_no_f = 1 + tr_tg_num - 1;           %adding the trial tag number of the first trial of this set to get the right filename
    trial_no_f = sprintf('%03.0f', trial_no_f);
    file_path = [direc f_name trial_no_f, '.tif'];  %full file path for currently analysed trial
    
    
    %skipping if cells have already been counted for dataset
    if exist([save_path 'n_cells_counted.txt']) == 2
        continue
    else
    end
    
    
    %reading in 10 frames of trial 1 to obtain an averaged image
    %for manual cell-counting
    for frame_n = 5:85
        frame = double(imread(file_path, frame_n));
        
        if frame_n == 5
            ave_frame = frame;
        else
        end
        
        if frame_n > 5
            ave_frame = ave_frame + frame;
        else
        end
    
    end
    ave_frame = ave_frame./10;
    disp_im = ave_frame;
    
    %manually selecting cells for count
    more = 1;
    counter = 0;
    n_skip_to = 1;
    maxp = max(max(disp_im)).*0.8;
    minp = min(min(disp_im));
    while more == 1
        counter = counter + 1;
        
        %plotting averaged image and manually picking the centre of a cell
        %to mark out and count
        colormap('gray');
        
        figure(1)
        imagesc(disp_im, [minp, maxp])
        set(gcf, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);
        
        [x, y] = ginput(1);
        
        disp_im = draw_square(disp_im, x, y, 2, maxp);
        
        figure(1)
        imagesc(disp_im, [minp, maxp])
        set(gcf, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);
        
        if counter == n_skip_to
            button = questdlg('More cells?','Select more cells?','~25 more','1 more', 'done', 'done');
        
        else
        end
        
        %parsing GUI input
        if strcmp(button, '~25 more') == 1
            n_skip_to = counter + 25;
            button = 'a';
        elseif strcmp(button, '1 more') == 1
            n_skip_to = counter + 1;
            button = 'a';
        elseif strcmp(button, 'done')
            more = 0;
        else
        end
            
        
    end
    n_cells = counter;
    
    %saving
    save([save_path 'n_cells_counted.txt'], 'n_cells', '-ASCII')
    
end
fclose(fid);



