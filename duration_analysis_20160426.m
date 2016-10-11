clear all
close all

%list_direc = 'D:\Data\CSHL\dataset_list_stim_dur_20160316.txt'; %KCs all
%list_direc = 'D:\Data\CSHL\dataset_list_sustained_MB418B_20160414.txt'; %KC A'B'
%list_direc = 'D:\Data\CSHL\dataset_list_sustained_MB185B_20160426.txt'; %KC AB
%list_direc = 'D:\Data\CSHL\dataset_list_sustained_MB131B_20160515.txt'; %KC G
list_direc = 'D:\Data\CSHL\dataset_list_PN_GH146_20161002.txt'; %PN axons


color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');
a = colormap('bone');
greymap = flipud(a);
colormap(greymap)


suppress_plots = 1;       %1 - doesn't plot stuff, 0 - plots stuff
fid = fopen(list_direc);
direc_counter = 0;

saved_cell_data = cell(1, 45);

%figure property initialisation variables
plot_height = 200;
plot_width = 200;
axis_font_size = 8;
tick_length = [0.02, 0.05];
line_widths = 2;
marker_sizes = 1;
marker_sizes_f = 4;
box_l_width = 0.5;


%loop to go through all experiment datasets listed in list file
while 1
    direc_counter = direc_counter + 1;
    direc = fgetl(fid);
    
    if ischar(direc) ~= 1
        break
    else
    end
    
    %skipping a particular OK107 dataset where n_frames was < stim time
    del = strfind(direc, '20160310');
    if isempty(del) == 0
        continue
    else
    end
    
    
    direc = [direc, '\'];
    
    %replacing C: with D:
    a = strfind(direc, 'C:');
    direc(a) = 'D';
    
    %skipping if this dataset has already been analysed
    if exist([direc 'cell_classifications.mat']) == 2
        disp('dataset already analysed; skipping...')
        %continue
    else
    end
    
    
    %loading extracted raw fluorescence data matrices written by
    %raw_dff_extractor
    dataset = load([direc 'expt.mat']);
    dataset = dataset.data;
    raw_data_mat = load([direc 'expt_raw_traces.mat']);
    raw_data_mat = raw_data_mat.raw_data_mat;
    
    
    %calculating dF/F traces and creating the sparse, 4-D, nan-filled
    %dff_data_mat 
    [dff_data_mat, stim_mat, prot_switch_trials] = cal_dff_traces_20160317(raw_data_mat, dataset, list_direc);
    clear raw_data_mat
    odor_list = unique(stim_mat(:, 1) );
    del = find(odor_list == 0);
    odor_list(del) = [];
    n_odors = length(odor_list);
    n_frames = size(dff_data_mat, 1);
    n_cells = size(dff_data_mat, 2);
    n_trials = size(dff_data_mat, 3);
    stim_time = dataset(1).stim.stimLatency + 0.54;     %odor stim delay measured with PID is 540 ms
    frame_time = dataset(1).info.framePeriod;
    stim_frame = floor(stim_time./frame_time);
    an_trial_window = nan;
    odor_dur_list = unique(stim_mat(:, 2) );
    n_reps = dataset(1).stim.reps;
    
    %identifying sig responses on a single trial basis, and then sig
    %responder cells in any individual block
    [resp_areas, sig_trace_mat, sig_cell_mat, sig_cell_block_mati] = ...
        cal_sig_responses_20160404(dataset, dff_data_mat, stim_mat, prot_switch_trials, list_direc, an_trial_window);
        
    del = isnan(resp_areas(1, :));
    bad_tr_list = find(del == 1);                   %list of trials thrown away due to movement
    good_tr_list = 1:n_trials;
    good_tr_list(bad_tr_list) = [];
    n_prot_blocks = length(prot_switch_trials);     %number of protocol blocks in this dataset
    
    %making sure that at least half the trials in each block, for each odor were acquired
    curr_window = n_trials./(length(odor_list).*n_prot_blocks);
    bad_od_blk_counter = 0;
    for odor_n = 1:length(odor_list)
        odor_ni = odor_list(odor_n);
        for odor_dur_n = 1:length(odor_dur_list)
            odor_dur_ni = odor_dur_list(odor_dur_n);
            for block_n = 1:n_prot_blocks
                od_blk_window_trs = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, odor_ni, odor_dur_ni, block_n, an_trial_window, 1);

                n_bad_trs = length(intersect(od_blk_window_trs, bad_tr_list));            %number of bad trials in current block
                if n_bad_trs > curr_window/2
                    disp('More than half the trials missing for one odor for one block.')
                    bad_od_blk_counter = bad_od_blk_counter + 1;
                else
                end

            end
        end
    end
    if bad_od_blk_counter > 1
        
        %continue
    else
    end
    
    
    %plotting ave response traces of all cells
    for odor_n = 1:length(odor_list)
       odor_ni = odor_list(odor_n);
       for odor_dur_n = 1:length(odor_dur_list)
           odor_dur_ni = odor_dur_list(odor_dur_n);
           try
               curr_color = color_vec(odor_ni, :);
           catch
               keyboard
           end
           stim_end_frm = ceil(stim_frame + odor_dur_ni./frame_time);
           %----------------- loops to walk through each odor, delivered for each duration
           
           tr_list = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, odor_ni, odor_dur_ni, block_n, an_trial_window, 1);               %list of trials of this odor, given for this duration
           ave_dff_mat = squeeze(dff_data_mat(:, :, tr_list, odor_ni) );
           ave_dff_mat = nanmean(ave_dff_mat, 3);          
           
           if suppress_plots == 0
               figure(1)
               imagesc(squeeze(ave_dff_mat)', [0, 4])
               xlabel('frame number')
               ylabel('cell number')
               set(gcf, 'Color', 'w')
               title(['fly number ' int2str(direc_counter) ', odor number ' int2str(odor_ni) ', odor duration' int2str(odor_dur_list(odor_dur_n))]);
               add_stim_shading(1, [stim_frame, stim_end_frm], 0.20, curr_color)
               set_xlabels_time(1, frame_time, 0.5)
               
               del = input('press enter');
           else
           end
       end
    end
   
    clear big_matrix
    %loop to analyse each cell across odors and durations
    for cell_n = 1:n_cells
        clear cell_data;
                
               
        for odor_n = 1:n_odors
            odor_ni = odor_list(odor_n);
            
            
            stim_fr_saved = [];
            stim_end_fr_saved = [];
            long_trace = [];
            prev_trace_length = 0;
            for odor_dur_n = 1:length(odor_dur_list);
                
                %initialising variable to store each trace for this cell
                if odor_n == 1 && odor_dur_n == 1
                    cell_data.traces = zeros(n_frames, n_odors, length(odor_dur_list), n_reps) + nan;
                else
                end
                
                odor_dur_ni = odor_dur_list(odor_dur_n);                        %actual odor duration in s
                stim_end_fr = ceil(stim_frame + (odor_dur_ni./frame_time) );    %frame number when stim ends
                tr_list = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, odor_ni, odor_dur_ni, block_n, an_trial_window, 1);
                                
                cell_data.traces(:, odor_n, odor_dur_n, 1:length(tr_list)) = squeeze(dff_data_mat(:, cell_n, tr_list, odor_ni) );
                
                stim_fr_saved = [stim_fr_saved; stim_frame];
                stim_end_fr_saved = [stim_end_fr_saved; stim_end_fr];
                
                curr_trace = cell_data.traces(:, odor_n, odor_dur_n);
                
                stim_frs_saved(odor_dur_n, :) = [10, ((stim_end_fr - stim_frame) + 10)] + prev_trace_length;
                
            end
            
            cell_data.stim_start_frs = stim_fr_saved;
            cell_data.stim_end_frs = stim_end_fr_saved;
            
            
                      
            
        end
        
        %calling cell-classifier
        disp(['classifying cell ' int2str(cell_n) 'of ' int2str(n_cells)])
        [r_vecs, h_vecs] = cell_classifier_20160514(cell_data, frame_time);
        
        %logging classification results
        classifier_r_saved(:, cell_n, 1:n_odors) = r_vecs';
        classifier_h_saved(:, cell_n, 1:n_odors) = h_vecs';
        
    end
    
    %saving data to file for current dataset
    classific_data.r_vecs = classifier_r_saved;
    classific_data.h_vecs = classifier_h_saved;
    
    save([direc 'cell_classifications.mat'], 'classific_data');
    
    disp(['Analysing direc ' int2str(direc_counter + 1)])
    
end
fclose(fid);

