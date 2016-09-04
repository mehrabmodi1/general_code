function [] = cascade_plot(data_mat, ydim, offset, linew)
%cascade_plot syntax: cascade_plot(data_mat, ydim, offset, linew) cascade_plot plots 
%multiple traces in one window with a y-displacement for each row. 
%Each row is individually normalised to its maximum value. cascade_plot can 
%only recieve a 2-D matrix as input, and will consider 
%ydim as the y axis for the cascade plot. 'offset' is a multiplicative
%factor that changes the y-offset between traces. 'linew' sets the
%thickness of the lines used to make the plot.


% setting the default values for parameters in case of absence of input
if nargin == 1
    ydim = 1;
    offset = 1;
    linew = 1;
elseif nargin == 2
    offset = 1;
    linew = 1;
elseif nargin == 3
    linew = 1;
end

if length(size(data_mat)) > 2
    error('Data_mat has > 2 dimensions. Only 2-D matrix input to cascade_plot.')
else
end

if ydim == 2
    data_mat = data_mat';
else
end


for row_no = 1:size(data_mat, 1)
    max_val = max(data_mat(row_no, :));
    data_mat(row_no, :) = data_mat(row_no, :)./max_val + row_no.*1.2.*offset;
end

plot(data_mat', 'black', 'LineWidth', linew)




    
    
