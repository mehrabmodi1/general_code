clear all
close all

clear all
close all
list_direc = ['C:\Data\CSHL\dataset_list_glc_ara_20150917.txt']; %expt datasets
%save_direc = 'C:\Stuff\CSHL\Glenn lab\Paper_drafts\2015_Toshi_MBON_plasticity\Data_quality_followup\movie_frames';
%list_direc = ['C:\Data\CSHL\dataset_KC_DN_pairing_control_20150828.txt']; %control datasets  

play_movies = 1;            %toggling playing of movies for each dataset


if isempty(findstr(list_direc, 'KC_DN_pairing_control')) ~= 1
    set_type = 0;   %control set
else
    set_type = 1;   %expt set
end


fid = fopen(list_direc);
direc_counter = 0;


cell_list = ...
{
    [2, 3], [1, 4, 6];...
    [2, 3, 4, 5], [2];...
    [9, 2, 4], [2, 5, 9];...
    [2, 8, 13, 20], [1, 2, 3, 15];...
    [1, 8, 9, 10], [1, 2, 5];...
};



bad_cell_list_manual = ...
{
    [2], [4];...
    [], [];...
    [2], [2];...
    [], [2];...
    [9, 10], [1];...        %last two trials
};



%loop to go through all experiment datasets listed in list file
while 1
    direc_counter = direc_counter + 1;
    direc = fgetl(fid);
       
    if ischar(direc) ~= 1
        break
    else
    end
    
    direc = [direc '\'];
        
    %loading extracted raw fluorescence data matrices written by
    %raw_dff_extractor
    dataset = load([direc 'expt.mat']);
    dataset = dataset.data;
    raw_data_mat = load([direc 'expt_raw_traces.mat']);
    raw_data_mat = raw_data_mat.raw_data_mat;
    
    
    %calculating dF/F traces and creating the sparse, 4-D, nan-filled
    %dff_data_mat 
    [dff_data_mat, stim_mat, prot_switch_trials] = cal_dff_traces_20150710(raw_data_mat, dataset, list_direc);
    clear raw_data_mat
    
    n_frames = size(dff_data_mat, 1);
    n_cells = size(dff_data_mat, 2);
    n_trials = size(dff_data_mat, 3);
    stim_time = dataset(1).stim.stimLatency;
    frame_time = dataset(1).info.framePeriod;
    stim_frame = floor(stim_time./frame_time);
    

       
    %identifying sig responses on a single trial basis, and then sig
    %responder cells in any individual block
    [resp_areas, sig_trace_mat, sig_cell_mat, sig_cell_block_mati] = cal_sig_responses_20150825(dataset, dff_data_mat, stim_mat, prot_switch_trials, list_direc);
    
    tr_list = [3:10, 13:20];
    dropped_trials = isnan(resp_areas(1, :));
    dropped_trials = dropped_trials(tr_list);
    dropped_trials = find(dropped_trials == 1);
    tr_list(dropped_trials) = [];
    
    %building averaged frame movie
    movie_mat = zeros(255, 256, length(tr_list)) + nan;
    
    
    %checking if averaged frames have been previously extracted and saved to disk
    if exist([direc 'aved_tr_frames.mat']) ~= 2
        for tr_ni = 1:length(tr_list)
            trial_n = tr_list(tr_ni);
            curr_stack = dataset(trial_n).imageStack;
            baseline_frame = nanmean(curr_stack(:, :, (41 - 22):(41 - 2)), 3);      %41 is stim frame

            movie_mat(:, :, tr_ni) = baseline_frame;
            save([direc 'aved_tr_frames.mat'], 'movie_mat');
        end
    else
        movie_mat = load([direc 'aved_tr_frames.mat']);
        movie_mat = movie_mat.movie_mat;
        
    end
    
    odor_list = unique(stim_mat(:, 1));             %list of odor stimuli used in this dataset 
    ROI_mat_orig = dataset(1).ROI(2).roi;
    
    %identifying centroids of each ROI
    centroids = regionprops(ROI_mat_orig, 'Centroid');
    centroids = round(cat(1, centroids.Centroid));
            
    %going through each ROI and looking for changes in pixel profile across
    %trials.
    for ROI_n = 1:n_cells
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
        
        %creating a sampling neighbourhood radially outwards from ROI centroid
        curr_ROI_mat = zeros(size(ROI_mat_orig, 1), size(ROI_mat_orig, 2));
        curr_ROI_mat(curr_centroid(2), curr_centroid(1)) = 1;
        se = strel('disk', 8);
        curr_ROI_mat = imdilate(curr_ROI_mat, se);
        curr_pix = find(curr_ROI_mat == 1);
        
        clear se
        clear curr_ROI_mat
                
        %calculating corrcoef for each ROI's pixels across trials
        ROI_pix_mat = zeros(21, 21, length(tr_list));
        ROI_all_pix_mat = zeros(21, 21, length(tr_list));
        c_vec = zeros(1, length(tr_list));
        for trial_n = 1:length(tr_list)
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
        
        
        
        if min(c_vec) < 0.82
            keyboard
            %playMovie(ROI_all_pix_mat)
            
        else
        end
        
    end
    
    keyboard
    
    
end
fclose(fid);
    