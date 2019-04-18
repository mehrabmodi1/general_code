%Written by Mehrab N Modi, 20140915. Segments image into overlapping
%squares and clusters pixels based on correlations in Ca time series.


clear all
close all


%House-keeping
%-----------------------

direc = 'C:\Data\CSHL\20140908\odor_stim_set2\';                        %directory with raw data in it
save_direc = 'C:\Data\CSHL\Analysed_data\';                             %main directory to save processed data files

[no_frames, no_trials, frame_rate, zoom, d_path, f_name, tr_tg_no] = tiff_info(direc);
frame_time = 1./frame_rate;                                             %in s
cell_dia = 18;                                                         %in pixels (zoom 2.6)

save_path = [save_direc d_path];                                        %final sub-directory to save processed data files
mkdir(save_path);

%generating mean fluorescence image
trial_no = 4;
trial_no_f = trial_no + tr_tg_no - 1;           %adding the trial tag number of the first trial of this set to get the right filename
trial_no_f = sprintf('%03.0f', trial_no_f);
filepath = [direc f_name trial_no_f, '.tif'];  %full file path for currently analysed trial
for frame_no = 1:30
    frame = imread(filepath, frame_no);
    frame = double(frame);
    if frame_no == 1
        ave_frame = zeros(size(frame));
    elseif frame_no > 1
        ave_frame = ave_frame + frame;
    else
    end
    
end
ave_frame = ave_frame./30;
n_cols = size(ave_frame, 2);
n_rows = size(ave_frame, 1);

%constructing segment framework
cell_dia_multiplier = 2;
seg_length = cell_dia.*cell_dia_multiplier;
k = cell_dia_multiplier.^2;                       %setting clustering k value based on how many cells to expect per segment
n_segs_c = ceil(n_cols./seg_length);
n_segs_r = ceil(n_rows./seg_length);
c_iters = 1;                                    %number of iterations of kmeans++ clustering algorithm


raw_data_mat = zeros(n_rows, n_cols, no_frames) + nan;

%loop to read in multiple trials
%for trial_no = 1:3
    trial_no = 4;
    trial_no_f = trial_no + tr_tg_no - 1;           %adding the trial tag number of the first trial of this set to get the right filename
    trial_no_f = sprintf('%03.0f', trial_no_f);
    file_path = [direc f_name trial_no_f, '.tif'];  %full file path for currently analysed trial
    
    %loop to read in all frames of current trial
    for frame_no = 1:no_frames
        frame = double(imread(filepath, frame_no));
        raw_data_mat(:, :, frame_no) = frame;
    end
%end


%loop to go through each image segment and calculate correlation matrix
ROI_mat = zeros(size(ave_frame));
for seg_c_no = 1:n_segs_c
    for seg_r_no = 1:n_segs_r
        r_no = seg_r_no;                                     %segment row number
        c_no = seg_c_no;                                     %segment col number

        %calculating co-ordinates of corners of segment in pixel row and col numbers
        lower_ri = ( ( (r_no - 1).*seg_length) + 1);         %lower row coord value
        lower_ci = ( ( (c_no - 1).*seg_length) + 1);         %lower col coord value
        
        %accounting for smaller segments at the edges
        if seg_r_no == n_segs_r
            higher_ri = n_rows;                                                  %higher row coord value
        else
            higher_ri = r_no.*seg_length;                                        %higher row coord value
        end

        if seg_c_no == n_segs_c
            higher_ci = n_cols;                                                   %higher col coord value
        else
            higher_ci = c_no.*seg_length;                                         %higher col coord value
        end
        n_pix = (higher_ri - lower_ri + 1).*(higher_ci - lower_ci + 1);           %number of pixels in each segment
        
        seg_pix_data = raw_data_mat(lower_ri:higher_ri, lower_ci:higher_ci, :);   %time series raw data for each pixel in current segment
        ave_seg = mean(seg_pix_data, 3);
        seg_pix_data = squeeze(reshape(seg_pix_data, n_pix, no_frames, 1));
        
        corr_mat = corrcoef(seg_pix_data');                                       %calculating corrcoef matrix
                
        
        %Clustering corrcoefs using kmeans_pp
        for c_iter = 1:c_iters
            [L, C, sum_dists] = kmeans(corr_mat, k, 'emptyaction', 'drop');
            
        end    
            %re-ordering corr_mat as per clustering identities
            pix_list = (1:1:size(corr_mat, 1))';
            pix_list = [L, pix_list];
            pix_list = sortrows(pix_list);
            
            corr_mat_ord = zeros(size(corr_mat));
            seg_pix_ord = zeros(size(seg_pix_data));
            
            for pix_n = 1:n_pix
                curr_pix = pix_list(pix_n, 2);
                corr_mat_ord(pix_n, :) = corr_mat(curr_pix, :);
                seg_pix_ord(pix_n, :) = seg_pix_data(curr_pix, :);
            end
 
            
            %visualising pixels according to cluster identity
            seg_ROI_mat = zeros(size(ave_seg));
            for clust_no = 1:k
                curr_pix = find(pix_list(:, 1) == clust_no);                    %indices of pixels in current cluster
                seg_ROI_mat(curr_pix) = clust_no;
                
                
                [curr_pix_r, curr_pix_c] = ind2sub(size(ave_seg), curr_pix);    %row-col subscripts of pixels in current cluster
                
                curr_pix_r = curr_pix_r + lower_ri;                             %adding offset for segment location within image frame
                curr_pix_c = curr_pix_c + lower_ci;                             %adding offset for segment location within image frame
                
                ROI_mat(curr_pix_r, curr_pix_c) = clust_no;
            end
            
            figure(1)
            subplot (2, 1, 1)
            imagesc(ave_seg)
            subplot(2, 1, 2)
            imagesc(seg_ROI_mat)
            
    end
end

figure(1)
subplot (2, 1, 1)
imagesc(ave_frame)
subplot(2, 1, 2)
imagesc(ROI_mat)

 