clear all
close all

list_direc = ['D:\Data\CSHL\dataset_list_sustained_MB418B_20160414.txt']; %expt datasets

choose_cell_type = 0;       %switch to turn on manual determination of cell-types.
type_to_plot = 2;           %select cell type to plot, only works if cell types have already been chosen.
                            %[1-sustained, 2-on, 3-off, 4-misc, 0-unreliable/nonresponders]
show_each_cell = 1;         %1-plot traces for interesting cells individually, 0 - skip to population plots
fid = fopen(list_direc);
direc_counter = 0;
sustained_counts = zeros(1, 2);
on_counts = zeros(1, 2);
off_counts = zeros(1, 2);
both_res_type_list_saved = [];

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
    
    direc = [direc '\'];
        
    %replacing C: with D:
    a = findstr(direc, 'C:');
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
                if n_bad_trs > curr_window./2
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
    
    
    %visualising response traces of each cell
    for odor_n = 1:length(odor_list)
       odor_ni = odor_list(odor_n);
       for odor_dur_n = 1:length(odor_dur_list)
           odor_dur_ni = odor_dur_list(odor_dur_n);           
           %----------------- loops to walk through each odor, delivered for each duration
           
           tr_list = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, odor_ni, odor_dur_ni, block_n, an_trial_window, 1);               %list of trials of this odor, given for this duration
           ave_dff_mat = squeeze(dff_data_mat(:, :, tr_list, odor_ni) );
           ave_dff_mat = nanmean(ave_dff_mat, 3);
%          
           


%            ave_dff_mat = zeros(n_frames, n_cells);
%            for tr_n = 1:length(tr_list)
%                tr_ni = tr_list(tr_n);
%                curr_dff_mat = dff_data_mat(:, :, tr_ni, odor_ni);
%                ave_dff_mat(:, :, 2) = curr_dff_mat;
%                ave_dff_mat = nansum(ave_dff_mat, 3);
%            end
%            ave_dff_mat = ave_dff_mat./length(tr_list);
           figure(1)
           imagesc(squeeze(ave_dff_mat)', [0, 4])
           xlabel('frame number')
           ylabel('cell number')
           set(gcf, 'Color', 'w')
           a = colormap('bone');
           a = flipud(a);
           colormap(a)
           title(['fly number ' int2str(direc_counter) ', odor number ' int2str(odor_ni) ', odor duration' int2str(odor_dur_list(odor_dur_n))]);
           
           
           %First make cell lists and then enable the lines for plotting below
           %cell_list = cell_list_maker(direc, odor_ni);      %done choosing cells manually
            
           del = input('press enter');
       end
    end 
       %Plotting
%        cell_type_list = [];
%        cell_list = load([direc 'cell_list_odor_' int2str(odor_ni) '.txt']);      %manually chosen interesting cells for this odor
% 
%        if show_each_cell == 1
%            for cell_n = 1:length(cell_list)
%                if choose_cell_type == 0
%                    cell_type_list = load([direc 'cell_type_list_odor_' int2str(odor_ni) '.txt']);
%                    curr_type = cell_type_list(cell_n);
%                    if curr_type == type_to_plot;
%                    else
%                        continue
%                    end
% 
%                else
%                end    
% 
% 
%                close(1)
%                stim_fr_saved = [];
%                stim_end_fr_saved = [];
%                for odor_dur_n = 1:length(odor_dur_list)
%                    odor_dur_ni = odor_dur_list(odor_dur_n);
%                    tr_list = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, odor_ni, odor_dur_ni, block_n, an_trial_window, 1);               %list of trials of this odor, given for this duration
%                    stim_end_time = stim_time + odor_dur_ni;
%                    stim_end_frame = ceil(stim_end_time./frame_time);
% 
%                    figure(1)
%                    fig_n = 1;
%                    subplot_vec = [1, 3, odor_dur_n];
%                    normf = 1;
%                    color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');
%                    stim_color = color_vec(odor_ni, :);               
% 
%                    title_st = ['fly ' int2str(direc_counter) ', odor ' int2str(odor_ni) ', duration ' int2str(odor_dur_ni) ', cell ' int2str(cell_n)];
%                    title_st = [];
%                    plot_traces_20160320(dff_data_mat, tr_list, stim_frame, stim_end_frame, stim_color, cell_n, odor_ni, fig_n, subplot_vec, normf, title_st)
%                    set_xlabels_time(1, frame_time, 0.3)
%                    hold on
%                    
%                    %saving traces to a mat file
%                    cell_.odor = odor_n;
%                    cell_.type = curr_type;
%                    
%                    stim_fr_saved = [stim_fr_saved; stim_frame];
%                    stim_end_fr_saved = [stim_end_fr_saved; stim_end_frame];
%                    
%                    if odor_dur_n == 1
%                        cell_.odor_dur1_traces = dff_data_mat(:, cell_n, tr_list, odor_ni);
%                    elseif odor_dur_n == 2
%                        cell_.odor_dur2_traces = dff_data_mat(:, cell_n, tr_list, odor_ni);
%                    elseif odor_dur_n == 3
%                        cell_.odor_dur3_traces = dff_data_mat(:, cell_n, tr_list, odor_ni);
%                    else
%                    end
%                    
%                end
%                cell_.stim_on_fr = stim_frame;
%                cell_.stim_off_fr = stim_end_frame;
%                
%                
%                    
%                fig_handle = figure(1);
%                set(fig_handle, 'Position', [25, 100, 1000, 225]);
%                del = input('cell type? sustained - 1, on - 2, off - 3, misc - 4, nothing - 0');
%                cell_type_list = [cell_type_list; del];
%            
% 
%            end
%                      
%        else
%        end
%        
%        if choose_cell_type ~=0
%            save([direc 'cell_type_list_odor_' int2str(odor_ni) '.txt'], 'cell_type_list', '-ASCII');
%        else
%        end
%        cell_type_list = load([direc 'cell_type_list_odor_' int2str(odor_ni) '.txt']);
%        sustained_count = length(find(cell_type_list == 1));
%        on_count = length(find(cell_type_list == 2));
%        off_count = length(find(cell_type_list == 3));
%        sustained_counts(1, odor_n) = sustained_counts(1, odor_n) + sustained_count;
%        on_counts(1, odor_n) = on_counts(1, odor_n) + on_count;
%        off_counts(1, odor_n) = off_counts(1, odor_n) + off_count;
%        
%        
%        %looking at cell-types across odors
%        if odor_n == 1
%            cell_list_saved = cell_list;
%            cell_types_saved = cell_type_list;
%        elseif odor_n == 2
%            [both_responders, cell_n_list2, cell_n_list1] = intersect(cell_list, cell_list_saved);
%            both_res_type_list = zeros(length(both_responders), 2);
%            both_res_type_list(:, 1) = cell_types_saved(cell_n_list1);
%            both_res_type_list(:, 2) = cell_type_list(cell_n_list2);
%            
%            %getting rid of cells that were marked manually as
%            %non-reliable responders
%            del = find(both_res_type_list(:, 1) == 0);
%            both_res_type_list(del, :) = [];
%            del = find(both_res_type_list(:, 2) == 0);
%            both_res_type_list(del, :) = [];
%           
%            both_res_type_list_saved = [both_res_type_list_saved; both_res_type_list];
%            
%        end
%        
%     end
%     
   
end
fclose(fid);

% figure(1)
% both_res_type_list_saved = sortrows(both_res_type_list_saved);
% imagesc(both_res_type_list_saved)
% ax = gca;
% ax.XTick = [1, 2];
% ax.XTickLabel = {'Odor3','Odor4'};
% ylabel('cell number')
% hcb = colorbar;
% hcb.YTick = [1:4];
% hcb.YTickLabel = {'sustained', 'on', 'off', 'mixed/ramp'};
% 
% 
% all_counts = [sustained_counts; on_counts; off_counts];
% figure(2)
% bar(all_counts)
% ylabel('number of cells')
% ax = gca;
% ax.XTickLabel = {'Sustained','On','Off'};
% 
% 
% 
% 
% 
% 


