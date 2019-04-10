function [BW_out] = stretch_maj_axis(BW, stretch_pix)
%This function takes a binary matrix with multiple objects, identifies the
%major axis of each object and stretches the object along its major axis by
%stretch_pix pixels.

if nargin < 2
    stretch_pix = 2;
else
end

% stretch_pix = 3;
% 
% BW1 = zeros(100, 100);
% BW1 = draw_circle(BW1, 50, 50, 5);
% line = strel('line', 8, 90);
% BW1 = imdilate(BW1, line);
% 
% BW2 = zeros(100, 100);
% BW2 = draw_circle(BW2, 25, 25, 5);
% line = strel('line', 8, 45);
% BW2 = imdilate(BW2, line);
% 
% BW = BW1 + BW2;
% 

labels = bwlabel(BW);
BW_out = zeros(size(BW));

for patch_n = 1:max(max(labels))
    curr_BW = zeros(size(BW));
    curr_patchi = find(labels == patch_n);
    curr_BW(curr_patchi) = 1;
    curr_BW = bwlabel(curr_BW);
    
    orien = regionprops(curr_BW, 'orientation');        %orientation of the current patch's maj axis
    orien = orien.Orientation;
    
    %stretching current patch along maj axis orientation
    line = strel('line', stretch_pix, orien);
    curr_BW = imdilate(curr_BW, line);
    
    BW_out = BW_out + curr_BW;
        
end


end