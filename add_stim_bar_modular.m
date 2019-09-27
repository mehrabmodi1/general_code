 function [] = add_stim_bar_modular(fig_n, stim_frames, stim_colors)
%Syntax: [] = add_stim_bar(fig_n, stim_frames, stim_colors)
%this function creates a fake, second pair of axes that is invisible, under
%the first one, and creates a patch colored stim_color in it to denote the stimulus duration.
%stim_colors is a matrix of color vecs with the number of rows equal to the
%number of patches.

%testing lines
% clear all
% close all
% test_vec = rand(1, 100);
% stim_frames = [20, 60];
% stim_colors = [.3, .5, .7];
% fig_n = 1;
% figure(1);
% plot(test_vec);

%checking if stim_frames has two sets of stim frames specified, one for each
%olfactometer or if it follows the old convention

if iscell(stim_frames) == 1
    if max(isnan(stim_frames{1, 2})) == 0
        n_olfs = 2;
    else
        n_olfs = 1;
    end
elseif iscell(stim_frames) == 0
    n_olfs = 1;
    stim_frames_orig = stim_frames;
    clear stim_frames
    stim_frames{1, 1} = stim_frames_orig;
end

disp('Run this function only AFTER fig_wrapup or the time-scale is distorted.');
figure(fig_n);
a = gca;
orig_ax = axis;
new_ax = [orig_ax(1:2), 0, 5];
orig_pos = a.Position;
new_pos = [orig_pos(1), (orig_pos(2) + orig_pos(4)), orig_pos(3), (orig_pos(4).*1.15) ];
ax2 = axes('Position', new_pos);    %created new axes object in same figure
axis(new_ax);
set(ax2, 'Color', 'none', 'Visible', 'off');
n_patches = size(stim_frames, 1);

%replicating stim_color and using it for all patches if only one color is
%specified for a pulse train
if n_olfs > 1 && size(stim_colors, 1) == 1
    stim_colors = repmat(stim_colors, n_olfs, 1);
else
end

for olf_n = 1:n_olfs
    curr_stim_frs = stim_frames{1, olf_n};
    offset = (olf_n - 1).*0.2;  %adding gap between rows of patches
    curr_colour = stim_colors(olf_n, :);
    
    for patch_n = 1:size(curr_stim_frs, 1)        %loop to draw multiple patches if needed
        %patch height
        y1 = .1 + offset;
        y2 = .2 + offset;

        %patch width
        x1 = curr_stim_frs(patch_n, 1);
        x2 = curr_stim_frs(patch_n, 2);

        %drawing patch
        y_vec = [y1, y2, y2, y1];
        x_vec = [x1, x1, x2, x2];

        p = patch(x_vec', y_vec', curr_colour );

        p.EdgeColor = [1, 1, 1];
        
    end
end