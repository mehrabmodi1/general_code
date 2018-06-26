function color_vec = setup_std_color_vec(n_colours)
%syntax: color_vec = setup_std_color_vec(n_colours)
%This function sets up a nice set of RGB triplets given the number of
%colours neeed.

color_vecs = load(['C:\Data\Code\general_code\std_color_vecs.mat']);
color_vecs = color_vecs.color_vecs;

if n_colours <= 6
    color_vec = color_vecs{1, 1};
elseif n_colors > 6 && n_colors <= 8
    color_vec = color_vecs{2, 1};
elseif n_colors > 8 %&& n_colors <=12
    color_vec = color_vecs{3, 1};
else
end
