function [reg_stack] = translate_stack(curr_stack, lags, pad)
%Syntax: [reg_stack] = translate_stack(curr_stack, lags, pad)
%This function translates all frames in a 3-D matrix in dims 1 and 2 by
%lags(1, 1) and lags(1, 2) places. Gaps created at the edges are padded with pad.

row_lag = round(lags(1, 1));
col_lag = round(lags(2, 1));
try 
    reg_stack = circshift(curr_stack, row_lag, 1);
    reg_stack = circshift(reg_stack, col_lag, 2);
catch
    keyboard
end

%replacing circularly shifted pixels with nans
if sign(col_lag) == 1
    reg_stack(:, 1:col_lag, :) = pad;
elseif sign(col_lag) == -1
    reg_stack(:, (end + col_lag):end, :) = pad;
else
end

if sign(row_lag) == 1
    reg_stack(1:row_lag, :, :) = pad;
elseif sign(row_lag) == -1
    reg_stack((end + row_lag):end, :, :) = pad;
else
end
