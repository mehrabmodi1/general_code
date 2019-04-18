
function [paper_colormap] = paper_colormap(length_vec)
%creates the custom, paper colormap


global paper_colormap


paper_colormap = [.3, .5, .85; .45, .85, .5; .80, .45, .4; .65, .9, .9];

paper_colormap = paper_colormap(1:length_vec, :);