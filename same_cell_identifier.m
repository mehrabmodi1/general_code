%This script automatically opens the ROI matrices for a pair of
%water-agarose datasets and allows the user to click on each image to
%identify the same cells across datasets. It then saves the paired-list of 'same
%cells' in both data directories.

clear all
close all

list_path = 'C:\Data\CSHL\dataset_list_20141031.txt';
save_direc = 'C:\Data\CSHL\Analysed_data\';

fid=fopen(list_path);
counter = 0;
while 1
    counter = counter + 1;
    direc = fgetl(fid);
    if ~ischar(direc),   break,   end
    direc = [direc '\'];
    
    
    disp(['dataset number ' int2str(counter)])
    
    if rem(counter, 2) == 1
        direc_w = direc;
    elseif rem(counter, 2) == 0
        direc_ag = direc;
        
        
        %------------------
        %reading in ROI matrices and ave_frames
        
        %reading in water dataset ROIs
        ROI_mat_w = load([direc_w 'ROI_mat.txt']);
        ROI_mat_W = im2bw(ROI_mat_w, 0.5);
        ROI_mat_w = bwlabel(ROI_mat_w);
        n_cells = max(max(ROI_mat_w));
        [n_frames, n_trials, frame_rate, zoom, d_path, f_name, tr_tg_num, ave_frame_w] = tiff_info(direc_w);
        save_path_w = [save_direc d_path];
        
        im_w = ave_frame_w;
        max_w = max(max(im_w));
        min_w = min(min(im_w));
        temp = find(ROI_mat_w > 0);
        im_w(temp) = max_w;
        
        
        %reading in agarose dataset ROIs
        ROI_mat_ag = load([direc_ag 'ROI_mat.txt']);
        ROI_mat_ag = im2bw(ROI_mat_ag, 0.5);
        ROI_mat_ag = bwlabel(ROI_mat_ag);
        n_cells = max(max(ROI_mat_ag));
        [n_frames, n_trials, frame_rate, zoom, d_path, f_name, tr_tg_num, ave_frame_ag] = tiff_info(direc_ag);
        save_path_ag = [save_direc d_path];
        
        im_ag = ave_frame_ag;
        max_ag = max(max(im_ag));
        min_ag = min(min(im_ag));
        temp = find(ROI_mat_ag > 0);
        im_ag(temp) = max_w;
        
        
        %plotting both images and asking user to manually pick the same
        %cell in the two images
        pair_list = [];
        more = 1;
        while more == 1
            colormap('gray')
            
            %plotting
            figure(1)
            subplot(2, 2, 1)
            imagesc(im_w, [min_w, max_w])

            subplot(2, 2, 2)
            imagesc(ave_frame_w)
            
            subplot(2, 2, 3)
            imagesc(im_ag, [min_w, max_w])
            
            subplot(2, 2, 4)
            imagesc(ave_frame_ag, [min_w, max_w])
            
            set(gcf, 'units', 'normalized', 'outerposition', [0, 0, 1, 1]);
            hold on

            
            %acquiring cell pair
            subplot(2, 2, 1)
            imagesc(im_w, [min_w, max_w])
            [x_w, y_w] = ginput(1);
            im_w = draw_square(im_w, x_w, y_w, 2, max_w./2);
            imagesc(im_w, [min_w, max_w])
            cell_n_w = ROI_mat_w(round(y_w), round(x_w) );         %cell number in the water ROI mat
            
            
            subplot(2, 2, 3)
            imagesc(im_ag, [min_w, max_w])
            [x_ag, y_ag] = ginput(1);
            im_ag = draw_square(im_ag, x_ag, y_ag, 2, max_w./2);
            imagesc(im_ag, [min_w, max_w])
            cell_n_ag = ROI_mat_ag(round(y_ag), round(x_ag) );     %cell number in the agarose ROI mat
            
            
            %asking user if they want to mark out more cell-pairs
            button = questdlg('More cell-pairs?','Select more cell-pairs?','yes','no', 'cancel');
            
            if strcmp(button, 'yes') == 1
                more = 1;
            elseif strcmp(button, 'no') == 1
                more = 0;
            elseif strcmp(button, 'cancel')
                more = 0;
            else
            end
            
            pair_list = [pair_list; cell_n_w, cell_n_ag];
        end
        
    
    
    else
    end

    if rem(counter, 2) == 1
        
    elseif rem(counter, 2) == 0

        if exist([save_path_w 'matched_cells_list.txt']) == 2
            overwrite = input('overwrite existing saved pairs? 1 - yes, 0 - no')

            if overwrite == 0
                continue
            else

            end

        else
        end

        save([save_path_w 'matched_cells_list.txt'], 'pair_list', '-ASCII')
        save([save_path_ag 'matched_cells_list.txt'], 'pair_list', '-ASCII')
    else
    end

end
fclose(fid);