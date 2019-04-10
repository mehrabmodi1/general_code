function [resp_areas, sig_trace_mat, sig_cell_mat] = cal_sig_responses(dataset, dff_data_mat, stim_mat)
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
n_odors = size(dff_data_mat, 4);

stim_time = dataset(1).stim.stimLatency.*1000;              %stimulus onset time in ms
stim_time = stim_time + 625;                                %added delay from valve opening to odor at pipe outlet
frame_time = dataset(1).info.framePeriod .* 1000;           %frame time in ms
stim_frame = floor(stim_time./frame_time);                  %frame no at which odor reached fly
pre_stim_frame = floor((stim_time - 4000)./frame_time);     %frame no 4 s prior to odor - to identify a 4s baseline

resp_areas = zeros(n_cells, n_trials) + nan;
sig_trace_mat = zeros(n_cells, n_trials);
sig_cell_mat = zeros(n_cells, n_odors) + nan;
for trial_n = 1:n_trials
    for cell_n = 1:n_cells
        curr_trace = squeeze(nanmean(dff_data_mat(:, cell_n, trial_n, :), 4));     %all traces but one along dim 4 are nan's, so it makes no difference
        base_trace = curr_trace(pre_stim_frame:(stim_frame - 2) );                 %pre-odor-stim baseline frames
        [del, pk_frame] = max(curr_trace(stim_frame:(stim_frame + ceil(1000./frame_time)) ));
        pk_frame = pk_frame + stim_frame;
        end_frame = pk_frame + round(3000./frame_time);
        if end_frame > n_frames
            end_frame = n_frames;
        else
        end
                
        resp_trace = curr_trace( (pk_frame - floor(300/frame_time)):end_frame ) ;
        resp_areas(cell_n, trial_n) = sum(resp_trace);
        
        base_m = mean(base_trace);
        resp_m = mean(resp_trace);
        base_s = std(base_trace);
        
        if resp_m > 4.*abs(base_m) && resp_m > 2.* base_s              %manually tweaked, p <<0.01 on a t-test with such a criterion
            h = 1;
            sig_trace_mat(cell_n, trial_n) = 1;
        else
            h = 0;
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

%identifying significantly responsive cells for each odor
for odor_n = 1:n_odors
    odor_trs = find(stim_mat(:, 1) == odor_n);                       %list of trials on which odor_n was delivered
    n_od_trs = length(odor_trs);
    for cell_n = 1:n_cells
        n_sig_resps = sum(sig_trace_mat(cell_n, odor_trs));    %no. of significant single trial responses of this cell to this odor
        
        if (n_sig_resps./n_od_trs) > 0.5                           %criterion as per Glenn's paper
            sig_cell_mat(cell_n, odor_n) = 1;
            
        else
            sig_cell_mat(cell_n, odor_n) = 0;
        end
    end
end



end