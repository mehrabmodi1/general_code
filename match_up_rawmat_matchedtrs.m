function [raw_data_mat_padded] = match_up_rawmat_matchedtrs(raw_data_mat, stim_mat_simple, tif_n_col_n)
%This function accounts for entries in stim_mat_simple (eg. for pairing
%trials) that have no corressponding, acquired .tiff but still need to be
%maintained in stim matrices. It inserts a dummy trial into raw_data_mat at
%the appropriate places.

matched_tif_ns = stim_mat_simple(:, tif_n_col_n);
nan_trs = find(isnan(matched_tif_ns) == 1);
raw_data_mat_padded = [];
pad = zeros(size(raw_data_mat, 1), size(raw_data_mat, 2), 1) + nan;
if isempty(nan_trs) == 0
    nan_tr_n1 = 0;      %just for the first nan_tr
    for nan_tr_n = 1:length(nan_trs)
        nan_tr_n2 = nan_trs(nan_tr_n);          %current NaN trial
        curr_matched_tif_list = matched_tif_ns( (nan_tr_n1 + 1):(nan_tr_n2 - 1) );     %list of matched .tif numbers from previous NaN trial to current NaN trial in stim_mat 
        raw_data_mat_padded = cat(3, raw_data_mat_padded, raw_data_mat(:, :, curr_matched_tif_list) );    %concatenating real data trials
        raw_data_mat_padded = cat(3, raw_data_mat_padded, pad);                                 %concatenating dummy trial
        nan_tr_n1 = nan_tr_n2;          %getting ready for next NaN trial, if any.
    end
    
    %concatenating remaining matched trials after last NaN trial
    curr_matched_tif_list = matched_tif_ns( (nan_tr_n2 + 1):length(matched_tif_ns) );
    raw_data_mat_padded = cat(3, raw_data_mat_padded, raw_data_mat(:, :, curr_matched_tif_list) );    %concatenating real data trials    
    
elseif isempty(nan_trs) == 1
    raw_data_mat_padded = raw_data_mat(:, :, stim_mat_simple(:, tif_n_col_n));     %case where every trial in stim_mat has a corressponding .tif acquired with it.
else
end