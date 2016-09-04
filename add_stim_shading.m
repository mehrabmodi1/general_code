function [] = add_stim_shading(figure_n, stim_frames, transparency, stim_color)
%This function adds a shaded, semi-transparent rectangular patch(es) over 
%figure_n, at x-axis intervals specified as pairs in stim_frames (each row 
%is a patch), with transparency and color as specified.

figure(figure_n)
hold on
axis_m = axis;
y1 = axis_m(1, 3);
y2 = axis_m(1, 4);

n_patches = size(stim_frames, 1);
%when only one color has been provided, all patches are made that color
if size(stim_color, 1) == 1
    stim_color = repmat(stim_color, n_patches, 1);
else
end

for patch_n = 1:n_patches
    curr_color = stim_color(patch_n, :);
    x1 = stim_frames(patch_n, 1);
    x2 = stim_frames(patch_n, 2);
    y_vec = [y1, y2, y2, y1];
    x_vec = [x1, x1, x2, x2];
    p = patch(x_vec', y_vec', curr_color);
    set(p, 'FaceAlpha', transparency);
    set(p, 'EdgeColor', 'none');
    set(gcf, 'Color', 'w')
    
end
axis(axis_m)
hold off