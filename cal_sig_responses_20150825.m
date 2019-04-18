function [resp_areas, sig_trace_mat, sig_cell_mat, sig_cell_block_mat] = cal_sig_responses(dataset, dff_data_mat, stim_mat, prot_switch_trials, list_direc)
%This function takes as inputs the sparse, 4-D dff_data_mat that contains dF/F traces stored
%by frame_n, cell_n, trial_n and odor_n; and stim_mat which contains
%stimulus delivery information for each trial. The outputs are three 2D
%matrices. Resp_areas of size n_cell, n_trials contains the area under the pk response 
%of each cell on each trial sig_trace_mat of size n_cells, n_trials indicates which
%individual traces were significant responses, and sig_cell_mat of size
%n_cells, n_odors which indicates which cells responded for more than half
%the presentations of any single odor.

n_frames = size(dff_data_mat, 1);
n_cells = size(dff_data_mat, 2);
n_trials = size(dff_data_mat, 3);
odor_list = unique(stim_mat(:, 1));
n_odors = length(odor_list);

stim_time = dataset(1).stim.stimLatency.*1000;              %stimulus onset time in ms
stim_time = stim_time + 625;                                %added delay from valve opening to odor at pipe outlet
frame_time = dataset(1).info.framePeriod .* 1000;           %frame time in ms
stim_frame = floor(stim_time./frame_time);                  %frame no at which odor reached fly
%pre_stim_frame = floor((stim_time - 4000)./frame_time);     %frame no 4 s prior to odor - to identify a 4s baseline
pre_stim_frame = 1;

resp_areas = zeros(n_cells, n_trials) + nan;
sig_trace_mat = zeros(n_cells, n_trials);
sig_cell_mat = zeros(n_cells, size(dff_data_mat, 4));
sig_cell_block_mat = zeros(n_cells, size(dff_data_mat, 4), length(prot_switch_trials));
for trial_n = 1:n_trials
    for cell_n = 1:n_cells
        curr_trace = squeeze(nanmean(dff_data_mat(:, cell_n, trial_n, :), 4));     %all traces but one along dim 4 are nan's, so it makes no difference
        base_trace = curr_trace(pre_stim_frame:(stim_frame - 2) );                 %pre-odor-stim baseline frames
       
        resp_trace = curr_trace(stim_frame:(stim_frame + ceil(5000/frame_time)) );
        resp_areas(cell_n, trial_n) = mean(resp_trace);
                
        base_m = mean(base_trace);
        resp_m = mean(resp_trace);
        base_s = std(base_trace);
        
        resp_trace_f = tsmovavg(resp_trace', 's', 5);                   %moving window average filtering
        
        if resp_m > 4.*abs(base_m) && resp_m > 2.* base_s              %manually tweaked, p <<0.01 on a t-test with such a criterion
            if max(resp_trace_f) > base_m + 2.33* base_s          %manually tweaked, p <<0.01 on a t-test with such a criterion
                h = 1;
                sig_trace_mat(cell_n, trial_n) = 1;
            else
                h = 0;
            end
        else
        end
%         if h == 1
%             figure(1)
%             plot(curr_trace)
%             title(int2str(resp_m./base_m))
%             keyboard
%         else
%         end
        
    end
end

%identifying significantly responsive cells for each odor, in any one block
%longer than 4 trials
for odor_n = 1:n_odors
    odor_ni = odor_list(odor_n);
    odor_trsi = find(stim_mat(:, 1) == odor_ni);                       %list of trials on which odor_n was delivered
        
    for block_n = 1:length(prot_switch_trials)
        [del, block_trs] = curr_block_ends(prot_switch_trials, block_n);
        
        
        
        %hack to analyse only the first four trials of block3 for Toshi's
        %paper
        if isempty(findstr(list_direc, 'Toshi')) == 0
            if block_n == 3
                block_trs = 13:20;
            else
            end
        end
                
        
        odor_trs = intersect(odor_trsi, block_trs);
        n_od_trs = length(odor_trs);
        
        %skipping of block is too short
        if length(odor_trs) < 4
            continue
        else
        end
        
        for cell_n = 1:n_cells
            n_sig_resps = sum(sig_trace_mat(cell_n, odor_trs));    %no. of significant single trial responses of this cell to this odor
            
             n_dropped_trs = length(find(isnan(resp_areas(1, odor_trs)) == 1));
             n_od_trs_nan = n_od_trs - n_dropped_trs;
            
            if (n_sig_resps./n_od_trs_nan) > 0.5                           %criterion as per Glenn's paper
                sig_cell_mat(cell_n, odor_ni) = 1;                         %a cell is called significantly responsive if it responds significantly to the given odor, in even one block of trials
                sig_cell_block_mat(cell_n, odor_ni, block_n) = 1;
                
            else
                
            end
        end
    end
    
    
end



