function [bad_cell_list, saved_lags] = bad_cell_catcher(dataset, direc, tr_list, supress)
%Syntax: [bad_cell_list] = bad_cell_catcher(dataset, direc)
%This function computes averaged frames from the baseline period for each
%trial. It then takes each cell-ROI's local neighbourhood and compares
%these pixels in each averaged frame with the first one by calculating
%corrcoefs. It then flags cells with even one mismatched ave frame for
%manual inspection. If supress == 1, the function doesn't play aligned
%frames for manual rejection of bad cells.
%Mehrab Modi, 20150929


    movie_mat = zeros(255, 256, length(tr_list)) + nan;
    stim_time = dataset(1).stim.stimLatency;
    frame_time = dataset(1).info.framePeriod;       %in s
    stim_frame = floor(stim_time./frame_time);
    
    
    %checking if averaged frames have been previously extracted and saved to disk
    if exist([direc 'aved_tr_frames.mat']) == 2
        %checking if movie_mat is out of date or if it has the same number
        %of frames as good trials just selected.
        movie_mat = load([direc 'aved_tr_frames.mat']);
        movie_mat = movie_mat.movie_mat;
        if size(movie_mat, 3) == length(tr_list)
            make_fresh_movie = 0;
        else
            make_fresh_movie = 1;
        end
    else
        make_fresh_movie = 1;
    end
    
    
    
    if make_fresh_movie == 1
        for tr_ni = 1:length(tr_list)
            trial_n = tr_list(tr_ni);
            curr_stack = dataset(trial_n).imageStack;
            baseline_frame = nanmean(curr_stack(:, :, (stim_frame - 22):(stim_frame - 2)), 3);     

            movie_mat(:, :, tr_ni) = baseline_frame;
            save([direc 'aved_tr_frames.mat'], 'movie_mat');
        end
    elseif make_fresh_movie == 0
        movie_mat = load([direc 'aved_tr_frames.mat']);
        movie_mat = movie_mat.movie_mat;
    else
        keyboard
    end
    
    ROI_mat_orig = dataset(1).ROI(2).roi;
    n_cells = max(max(ROI_mat_orig));
    n_trials = size(dataset, 2);
    saved_lags = zeros(n_cells, n_trials, 2) + nan;
    
    %identifying centroids of each ROI
    centroids = regionprops(ROI_mat_orig, 'Centroid');
    centroids = round(cat(1, centroids.Centroid));
            
    %going through each ROI and looking for changes in pixel profile across
    %trials.
    bad_cell_list = [];
    for ROI_n = 1:n_cells
        disp(['current cell is ' int2str(ROI_n) ' of ' int2str(n_cells) ' cells.']);
        curr_centroid = centroids(ROI_n, 1:2);
        
        % accounting for edge effects
        del = find(curr_centroid < 11);
        curr_centroid(del) = 11;         
        if curr_centroid(1) > (size(ROI_mat_orig, 1) - 10)
            curr_centroid(1) = size(ROI_mat_orig, 1) - 10;
        else
        end
        if curr_centroid(2) > (size(ROI_mat_orig, 2) - 10)
            curr_centroid(2) = size(ROI_mat_orig, 2) - 10;
        else
        end
                        
        %calculating corrcoef for each ROI's pixels across trials
        ROI_pix_mat = zeros(21, 21, length(tr_list));
        ROI_all_pix_mat = zeros(21, 21, length(tr_list));
        c_vec = zeros(1, length(tr_list));
        for trial_n = 1:length(tr_list)
            curr_frame = movie_mat(:, :, trial_n);        %averaged frame for current trial
            
            ROI_all_pix_mat (:, :, trial_n) = curr_frame((curr_centroid(2) - 10):(curr_centroid(2) + 10), (curr_centroid(1) - 10):(curr_centroid(1) + 10) );
            ROI_all_pix_mat(:, :, trial_n) = ROI_all_pix_mat(:, :, trial_n)./max(max(ROI_all_pix_mat(:, :, trial_n)));          %normalising 
            
            %calculating correlation coefficients of each trial with first
            %trial for local neighbourhood around ROI's centroid.
%             tr1 = ROI_pix_mat(:, :, 1);
%             c_vec(1, trial_n) = mat_corrcoef(tr1, ROI_pix_mat(:, :, trial_n));
            
        end
        
        [stack_al, c_vec, lags] = align_stack(ROI_all_pix_mat, 1);       %calling a function that re-aligns each cell's tile 
        saved_lags(ROI_n, tr_list, 1:2) = lags;
        
        ROI_all_pix_mat_orig = ROI_all_pix_mat;
        ROI_all_pix_mat = stack_al;
        

        if min(c_vec) < 0.82
            repeat = 1;
            while repeat == 1
                fig_handle = figure('name', ['cell number' int2str(ROI_n)]);
                set(fig_handle, 'Position', [100, 100, 300, 300]);
                playMovie(ROI_all_pix_mat, 0.05, 0)
                %keyboard
                repeat = input('1 - replay movie; 0 - stop playing');
            end
            close all
            keep_cell = input('Keep cell?; 1 - yes, 0 - no');
            if keep_cell == 0
                bad_cell_list = [bad_cell_list; ROI_n];
            else
            end
        else
        end
    
    end

end