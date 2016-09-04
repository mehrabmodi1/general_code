function [c lags] = xcorr_cust(vec1, vec2, max_lags)
%xcorr_cust syntax: [c lags] = xcorr_cust(vec1, vec2, max_lags)
%This function uses Matlab's corrcoef function to calculate normalised
%correlation coefficients, with lags between vec1 and vec2. If vec1 and
%vec2 are identical, the output vector c is the autocorrelation, if
%different, it is the cross correlation. The lags used are given as an
%output in lags - useful for plotting the auto or cross correlation.
%Max_lags is an input that imposes a limit to the lags used for the
%calculation to 0 +/- max_lags. ie the number of lags used will be
%2.*max_lags + 1.

% vec1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
% vec2 = [1, 3, 5, 7, 9, 11, 13, 15, 17, 19];
% max_lags = max([length(vec1), length(vec2)]) - 1;

%initialising
if nargin <2
    error('for autocorrelation, give same vector twice as vec1 and vec2')
elseif nargin < 3
    max_lags = max([length(vec1), length(vec2)]) - 1;
else
end


%house-keeping
size1 = size(vec1);
if size1(1, 1) > 1 && size1(1, 2) > 1
    error('Only vectors accepted as inputs for vec1 and vec2.')
elseif size1(1, 2) > 1
    vec1 = vec1';
else
end

size1 = size(vec2);
if size1(1, 1) > 1 && size1(1, 2) > 1
    error('Only vectors accepted as inputs for vec1 and vec2.')
elseif size1(1, 2) > 1
    vec2 = vec2';
else
end
clear size1


%building vector of lags to be used
lags = [-max_lags:1:0, 1:1:max_lags];
 
 
%loop to calculate correlation coefficients
c = zeros(max_lags, 1) + nan;
for lagi = 1:length(lags)
     lag = lags(lagi);
     pad = zeros(abs(lag), 1);
     if sign(lag) == -1
         veca = [vec1; pad];
         vecb = [pad; vec2];
     elseif sign(lag) == 1
         veca = [pad; vec1];
         vecb = [vec2; pad];
     else
     end
     
     c_val = corrcoef(veca, vecb);
     c(lagi, 1) = c_val(1, 2);
    
end
clear lag
clear c_val