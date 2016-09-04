function [od_blk_window_trs] = odor_trial_list_builder_20160317(stim_mat, prot_switch_trials, odor_ni, odor_dur_ni, block_n, window, front)
%this function's output is a list of trial numbers in the dataset that
%correspond to first or last (specified by front = 1/0) window trials in a
%protocol block for a given odor number.
%syntax: [od_blk_tr_list] = odor_trial_list_builder(stim_mat, prot_switch_trials, odor_ni, window, front)

odor_tr_lista = find(stim_mat(:, 1) == odor_ni);
odor_dur_tr_list = find(stim_mat(:, 2) == odor_dur_ni);         %list of all trials of odor_dur_ni duration
odor_tr_list = intersect(odor_tr_lista, odor_dur_tr_list);      %list of trials of odor_ni, of duration odor_dur_ni

block_ends = curr_block_ends(prot_switch_trials, block_n);
del = find( odor_tr_list >= block_ends(1) & odor_tr_list <= block_ends(2) );
block_odor_trs = odor_tr_list(del);


if isnan(window) == 1
    od_blk_window_trs = block_odor_trs;
else 
    if front == 1
        block_odor_trs = sort(block_odor_trs, 'ascend');
        try
            od_blk_window_trs = block_odor_trs(1:window);
        catch
            keyboard
        end
    elseif front == 0
        block_odor_trs = sort(block_odor_trs, 'descend');
        od_blk_window_trs = block_odor_trs(1:window);
        od_blk_window_trs = sort(od_blk_window_trs, 'ascend');
        
    else
        error('front must be 0 or 1.')
    end
    
end