function [diff_vecs] = get_dur_diffs(cell_data, dataset)
%This function calculates the difference in the area under the dF/F trace
%for each stimulus duration with the stimulus duration longer than it.
%Hence the ouput is a vector of length n-1 where n is the number of stim
%durations.


traces = cell_data.traces;
n_odors = size(traces, 2);
n_durs = size(traces, 3);
st_frs = cell_data.stim_start_frs;
end_frs = cell_data.stim_end_frs;
frame_time = dataset(1).info.framePeriod;
end_frs_off = end_frs + ceil(2.*(1./frame_time) );        %extending window of analysis beyond end of odor stim to include off responses

%making sure extended analysis window is not greater than the total number
%of frames available
del = find(end_frs > size(traces, 1) );
end_frs(del) = size(traces, 1);

dif_vecs = zeros(n_durs, n_odors);

for odor_n = 1:n_odors
    for dur_n = 2:n_durs
        trace1 = traces(st_frs(dur_n):end_frs(dur_n), odor_n, (dur_n - 1) );
        trace2 = traces(st_frs(dur_n):end_frs(dur_n), odor_n, (dur_n) );
        diff_vecs(dur_n, odor_n) = nansum(trace2 - trace1);
    end
end

end