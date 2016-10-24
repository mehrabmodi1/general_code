function [resp_areas, sig_trace_mat, sig_cell_mat, sig_cell_block_mat] = cal_sig_responses_20160404(dataset, dff_data_mat, stim_mat, prot_switch_trials, list_direc, an_trial_window)
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
odor_dur_list = unique(stim_mat(:, 2));
n_odor_durs = length(odor_dur_list);

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
    stim_duration = dataset(trial_n).stim.duration;
    if isnan(stim_duration) == 1
        continue
    else
    end
    stim_end_fr = ceil(stim_frame + ((stim_duration.*1000)./frame_time) );
    
    for cell_n = 1:n_cells
        curr_trace = squeeze(nanmean(dff_data_mat(:, cell_n, trial_n, :), 4));     %all traces but one along dim 4 are nan's, so it makes no difference
        base_trace = curr_trace(pre_stim_frame:(stim_frame - 2) );                 %pre-odor-stim baseline frames
                
        n_segments = ceil(stim_duration./5) + 1;        %number of 5s segments in stimulus response to analyse for sig responses; + 1 for off responses 
        seg_length = ceil(5000/frame_time);
        pk_times = zeros(1, n_segments) + nan;
        for segment_n = 1:n_segments
            if segment_n == 1
                seg_start_fr = stim_frame;
            else
                seg_start_fr = prev_seg_end_fr + 1;     %picking up where prev segment ended
            end
            
            seg_end_fr = seg_start_fr + seg_length;     %segment end determined by where segment starts
                        
            %making sure a given segment ends when the stimulus ends and
            %doesn't extend beyond it
            if sign(seg_start_fr - stim_end_fr) ~= sign(seg_end_fr - stim_end_fr)
                seg_end_fr = stim_end_fr;       %truncating segment so it ends at the end of the odor stimulus and doesn't extend beyond it
            
            else
            end
            
            %in case off period segment runs longer than n frames acquired
            if seg_end_fr > size(curr_trace, 1)
                seg_end_fr = size(curr_trace, 1);
            else
            end
            
            try
                resp_trace = curr_trace(seg_start_fr:seg_end_fr);
                resp_trace_f = tsmovavg_m(curr_trace', 's', 5);
                resp_trace_f = resp_trace_f(seg_start_fr:seg_end_fr);
            catch
                keyboard
            end
            
            %saving only on-response part of curve to response areas matrix
            if segment_n == 1
                resp_areas(cell_n, trial_n) = nanmean(resp_trace);
            else
            end

            %in case recorded trace is not long enough to allow analysing the extra segment
            if size(resp_trace, 1) < 5
                pk_times(segment_n) = [];
                continue
            else
            end
            
            base_m = mean(base_trace);
            base_s = std(base_trace);

            

    %       if resp_m > 4.*abs(base_m) && resp_m > 2.* base_s              
            if max(resp_trace_f) > base_m + 2.33* base_s          %p <0.05 on a t-test with such a criterion
                h = 1;
                sig_trace_mat(cell_n, trial_n) = 1;
                [del, pk_fr] = max(resp_trace_f);
                pk_fr = pk_fr + seg_start_fr - 1;
                pk_times(1, segment_n) = pk_fr.*frame_time;       %in ms

                %PICK UP THREAD HERE
                %fix problem with smoothing of short traces (1s responses)

                keyboard
                continue                                          %no need to analyse more segments if trace already found significant

            else
                h = 0;

            end

            if h == 1
                figure(1)
                plot(curr_trace)
                title(int2str(resp_m./base_m))
                keyboard
            else
            end
            prev_seg_end_fr = seg_end_fr;   %keeping track of segment end to determine seg_start_fr for next segment.
            keyboard
        end
    end
end

%identifying significantly responsive cells for each odor, in any one block
for odor_n = 1:n_odors
    for odor_dur_n = 1:n_odor_durs
        odor_ni = odor_list(odor_n);
        odor_dur_ni = odor_dur_list(odor_dur_n);

        for block_n = 1:length(prot_switch_trials)
            [block_od_trs] = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, odor_ni, odor_dur_ni, block_n, an_trial_window, 1);


            %skipping if block is too short
            if length(block_od_trs) < 2
                continue

            else
            end

            for cell_n = 1:n_cells
                n_sig_resps = sum(sig_trace_mat(cell_n, block_od_trs));    %no. of significant single trial responses of this cell to this odor

                n_dropped_trs = length(find(isnan(resp_areas(1, block_od_trs)) == 1));
                
                n_od_trs_nnan = length(block_od_trs) - n_dropped_trs;
                
                if (n_sig_resps./(length(block_od_trs) - n_od_trs_nnan) ) > 0.25                          %criterion as per Glenn's paper
                    sig_cell_mat(cell_n, odor_ni) = 1;                         %a cell is called significantly responsive if it responds significantly to the given odor, in even one block of trials
                    sig_cell_block_mat(cell_n, odor_ni, block_n) = 1;

                else

                end
            end
        end

    end
end



