function [area, pk_frame, p] = find_pk_area(vec, stim_frame, frame_time)
%This function takes a single dF/F response vector, and identifies the time
%of the peak response after the stimulus frame. It then sums the area 100
%ms before this frame and 400 ms after it. This gives the area under the
%dF/F curve at the time of peak response after stimulus delivery.
%Mehrab Modi, 20141015

if nargin < 3
    error('This function needs three inputs, the response vector, the stimulus frame number and the duration of each frame in ms.')
else
end

pre_pk_frame = round(200./frame_time);
post_pk_frame = round(600./frame_time);

[del pk_frame] = nanmax(vec(stim_frame:end));

sample = vec(max([stim_frame, (pk_frame - pre_pk_frame)]):min([length(vec), (pk_frame + post_pk_frame)]));            %sample around peak

area = nansum(sample)./length(sample);


%checking if response is significant
baseline_pts = vec( (stim_frame-16):(stim_frame - 1) );
response_pts = vec( max([stim_frame, (pk_frame - pre_pk_frame + stim_frame) ]):min([length(vec), (pk_frame + post_pk_frame + stim_frame) ]) );
[h, p] = ttest2(baseline_pts, response_pts);

end