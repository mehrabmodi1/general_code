clear all
close all

clear all
close all
list_direc = ['C:\Data\CSHL\dataset_list_Toshi_KC_DN_Led_20150708.txt']; %expt datasets
save_direc = 'C:\Stuff\CSHL\Glenn lab\Paper_drafts\2015_Toshi_MBON_plasticity\Data_quality_followup\movie_frames';
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
    
    for odor_ni = 1:length(odor_list)
        
        %identifying list of significant cells for this odor
        odor_n = odor_list(odor_ni);
        odor_trs = find(stim_mat(:, 1) == odor_n);
        odor_trs = intersect(tr_list, odor_trs);
        sig_cell_list = find(sig_cell_mat(:, odor_n) == 1);
        
        curr_bad_cells = cell_list{direc_counter, odor_ni};     %list of cells manually identified; numbered as per sig cell numbering
        
        cells_oi = sig_cell_list(curr_bad_cells);               %cell numbers of the manually identified cells as per whole dataset numbering
        
        %assembling matrix of ROIx color coded by number in curr_bad_cells
        ROIs_mat = zeros(size(movie_mat, 1), size(movie_mat, 2));
        for cell_oi_n = 1:length(cells_oi)
            curr_cell_n = cells_oi(cell_oi_n);                  %current cell number as per whole dataset cell numbering
            curr_cell_small_list_n = curr_bad_cells(cell_oi_n); %current cell number as per sig cell numbering
            
            curr_pix = find(ROI_mat_orig == curr_cell_n);       %pixels for ROI of current cell
            ROIs_mat(curr_pix) = curr_cell_small_list_n;        %labeling ROI as per cell number in sig cell numbering scheme for easy reference
        end
        
        odor_name = dataset(1).stim.odourNames(odor_n).odour;
        
        repeat = 1;
        
        ROI_list = unique(ROIs_mat);
        ROI_list(1) = [];               %removing the 0 listed as an ROI
        
        while repeat == 1
            for frame_n = 1:size(movie_mat, 3)
                curr_frame = movie_mat(:, :, frame_n);

                %re-scaling current frame to highest ROI value to match
                %colormaps
                curr_frame = curr_frame./max(max(curr_frame));
                curr_frame = curr_frame.*max(max(ROIs_mat));

                %replacing ROI pixels in curr_frame with ROI numbers
                for ROI_ni = 1:length(ROI_list)
                    ROI_n = ROI_list(ROI_ni);
                    curr_pix = find(ROIs_mat == ROI_n);
                    curr_frame(curr_pix) = ROI_n;
                end
                
                s_path = [save_direc 'fly' int2str(direc_counter) '\' odor_name ];
                mkdir(s_path)
                
                if play_movies == 1
                    h = figure(1)
                    imagesc(curr_frame)
                    title(int2str(frame_n))
                    colormap('gray')
                    colorbar('eastoutside')
                    
                    %pause(.05)
                    
                    print(h, '-dtiffn', [s_path '\frame_' int2str(frame_n) '.tif'])
                    
                else
                end

                
            end

            %repeat = input('0 to stop, 1 to repeat');
            repeat = 0;

        end
        
        
        
        
        
    end
    
    
    
    
end
fclose(fid);
    