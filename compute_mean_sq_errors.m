function [error_vec] = compute_mean_sq_errors(trace_mat, ref_vec)
%syntax:[error_vec] = compute_mean_sq_error(trace_mat, ref_vec)
%This function calculates the mean squared error of each vector in trace_mat
%with the reference vector rev_vec. The output is a vector of length = the
%number of traces. ie, the size of trace_mat in the dimension where
%size(ref_vec, dim) = 1.

[del, vec_dim] = max(size(ref_vec));

if vec_dim == 2
    trace_mat = trace_mat';
    ref_vec = ref_vec';
else
end

%loop to calculate MSQ for each vector in trace_mat
error_vec = zeros(1, size(trace_mat, 2)) + nan;
for trace_n = 1:size(trace_mat, 2)
    curr_trace = trace_mat(:, trace_n);
    %skipping nan traces
    if sum(curr_trace, 'omitnan') == 0
        continue
    else
    end
   
    error_vec(1, trace_n) = mean((ref_vec - curr_trace).^2, 'omitnan');
    
end


