function [im_matrix_out] = add_scale_bar(im_matrix_in, scale_micronsperpixel, barsize_microns, bar_colour)
%syntax: [im_matrix_out] = add_scale_bar(im_matrix_in, scale_micronsperpixel, barsize_microns, bar_colour)
%This function adds a scale bar to an image matrix. bar_colour can be 0 (bar
%will be matrix's min value) or 1 (bar will be matrix's max value).

%computing size of bar
bar_size = round((1./scale_micronsperpixel).*barsize_microns);

%computing location of bar
y1 = size(im_matrix_in, 2).*0.865;
y2 = size(im_matrix_in, 2).*0.9;

x2 = size(im_matrix_in, 1).*0.85;
x1 = x2 - bar_size;

if bar_colour == 1
    bar_val = max(max(im_matrix_in)).*1.05;
elseif bar_colour == 0
    bar_val = min(max(im_matrix_in)).*0.95;
else
end
    
im_matrix_out = im_matrix_in;
im_matrix_out(y1:y2, x1:x2) = bar_val;