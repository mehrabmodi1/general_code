function [stack_out, lag] = align_image_rows(stack_in, forced_lag)
%Syntax: [stack_out, lag] = align_image_rows(stack_in, forced_lag)
%This function corrects row lags due to imperfect turn-around time
%compensation in linear-scanned images. It computes a cross correlation
%between the odd and even rows and assigns a lag based on that.

stack_in_orig = stack_in;
stack_in = mean(stack_in, 3, 'omitnan');    %using stack mean to compute lags

odd_rows = stack_in(1:2:end, :);
odd_rows = reshape(odd_rows, 1, []);        %reshaped all odd rows into one long vector
eve_rows = stack_in(2:2:end, :);
eve_rows = reshape(eve_rows, 1, []);        %reshaped all even rows into one long vector

if isempty(forced_lag) == 1
    [r, lags] = xcorr(odd_rows, eve_rows, 40);

    [del, max_ri] = max(r);     %highest correlation value
    lag = lags(max_ri);         %corressponding lag
    clear del
else
    lag = forced_lag;
end

stack_out = zeros(size(stack_in_orig, 1), size(stack_in_orig, 2), size(stack_in_orig, 3)) + nan;

%inserting even rows un-shifted
eve_rows_all = stack_in_orig(2:2:end, :, :);
stack_out(2:2:end, :, :) = eve_rows_all;

%inserting even rows shifted by lag
odd_rows_shifted = stack_in_orig(1:2:end, 1:(size(stack_in, 2) - lag), :);
stack_out(1:2:end, (lag + 1):end, :) = odd_rows_shifted;
