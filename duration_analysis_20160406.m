close all

list_direc = 'D:\Data\CSHL\dataset_list_stim_dur_20160316.txt'; %expt datasets

choose_cell_type = 0;       %switch to turn on manual determination of cell-types.
type_to_plot = 4;           %select cell type to plot, only works if cell types have already been chosen.
                            %[1-sustained, 2-on, 3-off, 4-misc, 0-unreliable/nonresponders]
show_each_cell = 1;         %1-plot traces for interesting cells individually, 0 - skip to population plots
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

cell_counter = 0;
%loop to go through all experiment datasets listed in list file
while 1
    direc_counter = direc_counter + 1;
    direc = [fgetl(fid), '\'];
       
    if ischar(direc) ~= 1
        break
    else
    end
    
    %replacing C: with D:
    a = strfind(direc, 'C:');
    direc(a) = 'D';
    
    
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
    n_odors = length(odor_list);
    n_frames = size(dff_data_mat, 1);
    n_cells = size(dff_data_mat, 2);
    n_trials = size(dff_data_mat, 3);
    stim_time = dataset(1).stim.stimLatency + 0.54;     %odor stim delay measured with PID is 540 ms
    frame_time = dataset(1).info.framePeriod;
    stim_frame = floor(stim_time./frame_time);
    an_trial_window = nan;
    odor_dur_list = unique(stim_mat(:, 2) );
       
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
           %----------------- loops to walk through each odor, delivered for each duration
           
           tr_list = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, odor_ni, odor_dur_ni, block_n, an_trial_window, 1);               %list of trials of this odor, given for this duration
           ave_dff_mat = squeeze(dff_data_mat(:, :, tr_list, odor_ni) );
           ave_dff_mat = nanmean(ave_dff_mat, 3);          
           
           
           figure(1)
           imagesc(squeeze(ave_dff_mat)', [0, 4])
           xlabel('frame number')
           ylabel('cell number')
           set(gcf, 'Color', 'w')
           a = colormap('bone');
           a = flipud(a);
           colormap(a)
           title(['fly number ' int2str(direc_counter) ', odor number ' int2str(odor_ni) ', odor duration' int2str(odor_dur_list(odor_dur_n))]);

           del = input('press enter');
       end
    end
    
    type_to_plot = 1;
    
    %loop to analyse each cell across odors and durations
    
    for cell_n = 1:n_cells
        clear cell_data;
        
        cell_type_list = load([direc 'cell_type_list_odor_' int2str(odor_ni) '.txt']);
        if cell_n > length(cell_type_list)
            continue
        else
        end
        
        
        curr_type = cell_type_list(cell_n);
        
        %skipping non-responsive cells
%         if curr_type ~= type_to_plot
%             continue
%         else
%         end
        if curr_type == 0
            continue
        else
        end
        
        for odor_n = 1:n_odors
            odor_ni = odor_list(odor_n);
            
            
            stim_fr_saved = [];
            stim_end_fr_saved = [];
            for odor_dur_n = 1:length(odor_dur_list);
                odor_dur_ni = odor_dur_list(odor_dur_n);                        %actual odor duration in s
                stim_end_fr = ceil(stim_frame + (odor_dur_ni./frame_time) );    %frame number when stim ends
                                
                tr_list = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, odor_ni, odor_dur_ni, block_n, an_trial_window, 1);
                
                if curr_type == 1
                    cell_data.type = 'sustained';
                elseif curr_type == 2
                    cell_data.type = 'on';
                elseif curr_type == 3
                    cell_data.type = 'off';
                elseif curr_type == 4
                    cell_data.type = 'mixed';
                end
                
                cell_data.traces(:, odor_n, odor_dur_n) = nanmean(squeeze(dff_data_mat(:, cell_n, tr_list, odor_ni) ), 2);
                cell_data.sds(:, odor_n, odor_dur_n) = nanstd(squeeze(dff_data_mat(:, cell_n, tr_list, odor_ni) ), [], 2);
                
                
                stim_fr_saved = [stim_fr_saved; stim_frame];
                stim_end_fr_saved = [stim_end_fr_saved; stim_end_fr];
                
                
                
            end

            cell_data.stim_start_frs = stim_fr_saved;
            cell_data.stim_end_frs = stim_end_fr_saved;
            
        end
   
        %analysing current cell, curent odor
        %smoothing traces, moving window average
        traces = cell_data.traces;
        traces_s = filter([.2, .2, .2, .2, .2], 1, traces);

        %looking if resp trace changes with stim duration
        diff_vecs = get_dur_diffs(cell_data, dataset); 
        
        cell_data.type;
        
        
        
        %plotting duration differences against each other
        for odor_n = 1:n_odors
            figure(4)
            
            if curr_type == 1
                curr_color = 'r';
            elseif curr_type == 2
                curr_color = 'g';
            elseif curr_type == 3
                curr_color = 'b';
            else
                continue
            end
            sum_diffs = sum(diff_vecs(2:3, :));
            plot(sum_diffs(1), sum_diffs(2), ['O' curr_color])
            
            hold on
        end
        
        
        
        cell_counter = cell_counter + 1;
        saved_cell_data{1, cell_counter} = cell_data;
    end
    
   
end
fclose(fid);




