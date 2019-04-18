function [] = trajectory_plotter(data_mat, figure_n, save_direc, line_width, color, color_seq, pause_mode, marker_size1)
%This function takes a 3D dataset with a maximum length of 3 along dim2 and plots
%a 2 or 3D plot of the time-series data, assuming dim 1 is time. If
%save_direc is not empty, it saves frames for the addition of each point to the
%plot. If color_seq is 0, all traces are plotted in 'color' color. If
%color_seq is 1, each trace is a different shade of color, getting lighter
%trace by trace. If pause_mode is 1, the function waits for manual input
%before it plots the next trace.
%syntax: function [] = trajectory_plotter(data_mat, figure_n, save_direc, line_width, color, color_seq, pause_mode, marker_size)

n_frames = size(data_mat, 1);
n_cells = size(data_mat, 2);
n_trials = size(data_mat, 3);

col_dark_lim = color./3;
col_light_lim = color./max(color);

r_vec = col_dark_lim(1):( (col_light_lim(1) - col_dark_lim(1))./n_trials ):col_light_lim(1);
g_vec = col_dark_lim(2):( (col_light_lim(2) - col_dark_lim(2))./n_trials ):col_light_lim(2);
b_vec = col_dark_lim(3):( (col_light_lim(3) - col_dark_lim(3))./n_trials ):col_light_lim(3);

if color_seq == 1
    color_vec = [r_vec', g_vec', b_vec'];
elseif color_seq == 0
    color_vec = color;
else
end

for trial_n = 1:n_trials
    if length(marker_size1 > 1)
        
        marker_size = marker_size1(trial_n);
        
    elseif length(marker_size1) == 1
        marker_size = marker_size1;
    end
    col_n = rem(trial_n, size(color_vec, 1));
    if col_n == 0
        col_n = size(color_vec, 1);
    else
    end
    curr_color = color_vec(col_n, :);
    for frame_n = 1:n_frames
        figure(figure_n)
        if n_cells == 2                     %2-D plot implies 2 PCs
           plot(data_mat(frame_n, 1, trial_n), data_mat(frame_n, 2, trial_n), 'O', 'Color', curr_color, 'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'MarkerSize', marker_size) 
           if frame_n > 1
               %connecting last two points with a line
               if line_width > 0
                    line([data_mat( (frame_n-1), 1, trial_n); data_mat(frame_n, 1, trial_n)], [data_mat((frame_n-1), 2, trial_n); data_mat(frame_n, 2, trial_n)], 'Color', curr_color, 'LineWidth', line_width)
               else
               end
           
           else
           end
            
        elseif n_cells == 3                 %3-D plot
            try
                plot3(data_mat(frame_n, 1, trial_n), data_mat(frame_n, 2, trial_n), data_mat(frame_n, 3, trial_n), 'O', 'Color', curr_color, 'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color, 'MarkerSize', marker_size)
            catch
                keyboard
            end
            if frame_n > 1
                %connecting last two points with a line
                plot3([data_mat((frame_n-1), 1, trial_n); data_mat(frame_n, 1, trial_n)], [data_mat((frame_n-1), 2, trial_n); data_mat(frame_n, 2, trial_n)], [data_mat((frame_n-1), 3, trial_n); data_mat(frame_n, 3, trial_n)], 'Color', curr_color, 'LineWidth', line_width, 'MarkerFaceColor', curr_color, 'MarkerEdgeColor', curr_color)
            else
            end
            
        else
        end
        pause(0.05);
        if isempty(save_direc) == 0
            print([save_direc '\trial' int2str(trial_n) '_frame' int2str(frame_n)], 'jpg');
        else
        end
        hold on
    end
    if trial_n < n_trials
        if pause_mode == 1
            del = input('press any key for next trial');
        else
        end
    else
    end

end


end