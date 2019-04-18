clear all
close all

direc_lists_mat =  [{'D:\Data\CSHL\dataset_list_sustained_MB418B_20160414.txt'};... %KC A'B'
                    {'D:\Data\CSHL\dataset_list_sustained_MB185B_20160426.txt'};... %KC AB
                    {'D:\Data\CSHL\dataset_list_sustained_MB131B_20160515.txt'};... %KC G
                    {'D:\Data\CSHL\dataset_list_stim_dur_20160316.txt'};... %OK107 KCs
                    %{'D:\Data\CSHL\dataset_list_sustained_OK107xsytGCaMP6s_20161212.txt'};... %vert lobe, sytGCaMP
                    ]; 
                
dump_direcs = [{'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Analysis\data_dump\ApBp\'};...
                {'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Analysis\data_dump\AB\'};...
                {'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Analysis\data_dump\G\'};...
                {'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Analysis\data_dump\OK107\'}];
            
            
n_direc_lists = size(direc_lists_mat, 1);
                
                
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');

%a = colormap('bone');
%greymap = flipud(a);


%figure property initialisation variables
plot_height = 200;
plot_width = 280;
axis_font_size = 8;
tick_length = [0.02, 0.05];
line_widths = 2;
marker_sizes = 1;
marker_sizes_f = 4;
box_l_width = 0.5;

[del, odor_names] = xlsread('D:\Data\CSHL\odor_names_20161108.xls', 1);

%loop to go through all directory lists
for direc_list_n = 1:n_direc_lists

    list_direc = direc_lists_mat{direc_list_n, 1};
    dump_direc = dump_direcs{direc_list_n, 1};
    
    fid = fopen(list_direc);
    direc_counter = 0;
    saved_resp_vecs = [];
    saved_resp_vecs_1s = [];
    saved_resp_vec_ratio = [];
    all_sig_cells_mat = [];
    all_sig_cells_1s_mat = [];
    
    
    %parsing direc list path for name of direc list
    namei = findstr(list_direc, 'dataset_list');
    dir_list_name = (list_direc(namei:(end-4)));
        
    %loop to go through all experiment datasets listed in list file
    cell_counter = 0;
    while 1
        
        direc = fgetl(fid);

        if ischar(direc) ~= 1
            break
        else
        end
        
        direc_counter = direc_counter + 1;

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

        %loading extracted raw fluorescence data matrices written by
        %raw_dff_extractor
        dataset = load([direc 'expt.mat']);
        dataset = dataset.data;
        raw_data_mat = load([direc 'expt_raw_traces.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;


        %calculating dF/F traces and creating the sparse, 4-D, nan-filled
        %dff_data_mat 
        [dff_data_mat, stim_mat, prot_switch_trials] = cal_dff_traces_20160317(raw_data_mat, dataset, list_direc);
        
%         dump_direc = dump_direcs{direc_list_n};
%         save([dump_direc 'dFF_data_fly' int2str(direc_counter) '.mat'], 'dff_data_mat');
%         save([dump_direc 'stim_info_fly' int2str(direc_counter) '.mat'], 'stim_mat');
        
        clear raw_data_mat
        odor_list = unique(stim_mat(:, 1) );
        del = find(odor_list == 0);
        odor_list(del) = [];
        del = find(isnan(odor_list) == 1);
        odor_list(del) = [];
        n_odors = length(odor_list);
        n_frames = size(dff_data_mat, 1);
        n_cells = size(dff_data_mat, 2);
        n_trials = size(dff_data_mat, 3);
        prot_switch_trials = n_trials + 1;                  
        stim_time = dataset(1).stim.stimLatency + 0.54;     %odor stim delay measured with PID is 540 ms
        frame_time = dataset(1).info.framePeriod;
        stim_frame = floor(stim_time./frame_time);
        an_trial_window = nan;
        odor_dur_list = unique(stim_mat(:, 2) );
        del = find(odor_dur_list == 0);
        odor_dur_list(del) = [];
        del = find(isnan(odor_dur_list) == 1);
        odor_dur_list(del) = [];
        n_reps = dataset(1).stim.reps;
        stim_end_frames = [round(odor_dur_list./frame_time) + stim_frame];
        
        
        %identifying sig responses on a single trial basis, and then sig
        %responder cells in any individual block
        [resp_areas, sig_trace_mat, sig_cell_mat, sig_trace_mat_old, sig_cell_mat_old] = ...
            cal_sig_responses_20161024(dataset, dff_data_mat, stim_mat, prot_switch_trials, list_direc, an_trial_window);
            
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
        
        long_dur_n = find(odor_dur_list == 60);
        n_reps = ceil(n_trials./(n_odors.*length(odor_dur_list)));
        
        
        %initialising trace logging data structure
        for cell_count_n = 1:n_cells
            saved_data(cell_counter + cell_count_n).traces = zeros(n_frames, n_reps, n_odors, length(odor_dur_list)) + nan;
            saved_data(cell_counter + cell_count_n).info.sig_resps = zeros(n_odors, length(odor_dur_list)) + nan;
            saved_data(cell_counter + cell_count_n).info.stim_frs = zeros(2, length(odor_dur_list));
            saved_data(cell_counter + cell_count_n).info.frame_time = frame_time;
        end
        
        
        for odor_n = 1:length(odor_list)
            odor_ni = odor_list(odor_n);
            
            %identifying significant responders to 60s stimuli, for current odor
            sig_cells = find(sig_cell_mat(:, odor_ni, long_dur_n) == 1);
            od_trs = find(stim_mat(:, 1) == odor_ni);
            
            stim_end_fr = stim_frame + ceil(60./frame_time);
            an_end_fr = min([n_frames, ceil(stim_end_fr + 200./frame_time)]);
            
            %building list of trs of current odor for each dur
            for dur_n = 1:length(odor_dur_list)
                curr_dur = odor_dur_list(dur_n);
                stim_end_fr = stim_frame + ceil(odor_dur_list(dur_n)./frame_time);
                an_end_fr = min([n_frames, ceil(stim_end_fr + 200./frame_time)]);
                dur_trs = find(stim_mat(:, 2) == curr_dur);
                curr_trs = intersect(od_trs, dur_trs);
                
                for cell_n = 1:n_cells
                    %logging traces and info for each cell
                    saved_data(cell_counter + cell_n).traces(:, 1:length(curr_trs), odor_n, dur_n) = squeeze(dff_data_mat(:, cell_n, curr_trs, odor_ni));
                    saved_data(cell_counter + cell_n).info.sig_resps(odor_n, dur_n) = sig_cell_mat(cell_n, odor_ni, dur_n);
                    saved_data(cell_counter + cell_n).info.stim_frs(:, dur_n) = [stim_frame, stim_end_fr];
                    
                end
            end
            
        end
        cell_counter = cell_counter + n_cells;
        
    end

    %saving cell data to disk for the current direc list
    mkdir(dump_direc);
    save([dump_direc 'cell_data.mat'], 'saved_data');
    
    
end
