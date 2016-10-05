function [bad_cell_list] = ROI_z_movt_checker_20160505(dataset, direc)
%This function plays a zoomed-in movie of each cell, across trials for the
%user to decide if its stable or moving in z across the dataset, allowing
%the user to throw away cells that had z-movt.


%building averaged frame movie
movie_mat = zeros(255, 256, length(dataset)) + nan;
n_trials = length(dataset);


%checking if averaged frames have been previously extracted and saved to disk
if exist([direc 'aved_tr_frames.mat']) ~= 2
    for trial_n = 1:n_trials
        curr_stack = dataset(trial_n).imageStack;
        baseline_frame = nanmean(curr_stack(:, :, (41 - 22):(41 - 2)), 3);      %41 is stim frame

        movie_mat(:, :, trial_n) = baseline_frame;
        save([direc 'aved_tr_frames.mat'], 'movie_mat');
    end
else
    movie_mat = load([direc 'aved_tr_frames.mat']);
    movie_mat = movie_mat.movie_mat;

end

ROI_mat_orig = dataset(1).ROI(2).roi;
n_cells = max(max(ROI_mat_orig));

%identifying centroids of each ROI
centroids = regionprops(ROI_mat_orig, 'Centroid');
centroids = round(cat(1, centroids.Centroid));

%going through each ROI and looking for changes in pixel profile across
%trials.
bad_cell_list = [];
for ROI_n = 1:n_cells
    curr_centroid = centroids(ROI_n, 1:2);

    % accounting for edge effects
    curr_centroid(curr_centroid < 11) = 11;         
    if curr_centroid(1) > (size(ROI_mat_orig, 1) - 10)
        curr_centroid(1) = size(ROI_mat_orig, 1) - 10;
    else
    end
    if curr_centroid(2) > (size(ROI_mat_orig, 2) - 10)
        curr_centroid(2) = size(ROI_mat_orig, 2) - 10;
    else
    end

    %creating a sampling neighbourhood radially outwards from ROI centroid
    curr_ROI_mat = zeros(size(ROI_mat_orig, 1), size(ROI_mat_orig, 2));
    curr_ROI_mat(curr_centroid(2), curr_centroid(1)) = 1;
    se = strel('disk', 8);
    curr_ROI_mat = imdilate(curr_ROI_mat, se);
    curr_pix = find(curr_ROI_mat == 1);

    clear se
    clear curr_ROI_mat

    %calculating corrcoef for each ROI's pixels across trials
    ROI_pix_mat = zeros(21, 21, n_trials);
    ROI_all_pix_mat = zeros(21, 21, n_trials);
    c_vec = zeros(1, n_trials);
    for trial_n = 1:n_trials
        curr_frame = movie_mat(:, :, trial_n);        %averaged frame for current trial
        curr_im_pix = zeros(size(movie_mat, 1), size(movie_mat, 2));
        curr_im_pix(curr_pix) = curr_frame(curr_pix);
        ROI_pix_mat(:, :, trial_n) = curr_im_pix((curr_centroid(2) - 10):(curr_centroid(2) + 10), (curr_centroid(1) - 10):(curr_centroid(1) + 10) );
        ROI_pix_mat(:, :, trial_n) = ROI_pix_mat(:, :, trial_n)./max(max(ROI_pix_mat(:, :, trial_n)));          %normalising 
        ROI_all_pix_mat (:, :, trial_n) = curr_frame((curr_centroid(2) - 10):(curr_centroid(2) + 10), (curr_centroid(1) - 10):(curr_centroid(1) + 10) );
        ROI_all_pix_mat(:, :, trial_n) = ROI_all_pix_mat(:, :, trial_n)./max(max(ROI_all_pix_mat(:, :, trial_n)));          %normalising 

        %calculating correlation coefficients of each trial with first
        %trial for local neighbourhood around ROI's centroid.
        tr1 = ROI_pix_mat(:, :, 1);
        c_vec(1, trial_n) = mat_corrcoef(tr1, ROI_pix_mat(:, :, trial_n));

    end
    
    
    
    skip_set = 0;
    if min(c_vec) < 0.82
        play_again = 1;
        
        while play_again == 1
            figure(1)
            playMovie(ROI_all_pix_mat, 0.05, 0)
            
            qstring = 'What do you want to do?';
            btn_options = {'keep cell', 'discard cell', 'replay movie', 'skip dataset', 'save movie'};
            btn_n = bttnChoiseDialog(btn_options, 'Input Needed', 1, qstring);
            
            if btn_n == 1
                keep_cell = 1;
                play_again = 0;
            elseif btn_n == 2
                keep_cell = 0;
                play_again = 0;
            elseif btn_n == 3
                play_again = 1;
            elseif btn_n == 4
                keep_cell = 1;
                play_again = 0;
                skip_set = 1;
            elseif btn_n == 5
                play_again = 1;
                
                %saving movie
                save_direc = uigetdir('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\', 'Select folder to save movie in.');
                c = clock;
                c = fix(c);
                save_direc = [save_direc, '\ROI_movie_', int2str(c(4)), int2str(c(5)), int2str(c(6)), '.tif'];
                
                n_frames = size(ROI_all_pix_mat, 3);
                for frame_n = 1:n_frames
                    imwrite(ROI_all_pix_mat(:, :, frame_n),save_direc,'WriteMode','append')
                end
                                
            else
            end
            
        end

    else
        keep_cell = 1;
    end
    if keep_cell == 0
        bad_cell_list = [bad_cell_list; ROI_n];
    else
    end
    
    if skip_set == 1
        break
    else
    end
        
        
end
