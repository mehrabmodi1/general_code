color_vec_path = 'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt';

color_vec = load(color_vec_path);

n_colors = size(color_vec, 1);

r_vec = rand(1, 250);

for color_n = 1:n_colors
    r_vec = r_vec + 2;
    curr_color = color_vec(color_n, :);
    figure(1)
    plot(r_vec, 'LineWidth', 2, 'Color', curr_color)
    hold on
    edit_in = input('1 if OK, 0 to edit');
    
    while edit_in ~= 1
        disp(num2str(curr_color))
        new_color = input('enter new color as a 3 num vector');
        color_vec(color_n, :) = new_color;
        curr_color = new_color;
        figure(1)
        plot(r_vec, 'LineWidth', 2, 'Color', curr_color)
        edit_in = input('1 if OK, 0 to edit');
    end
    
end
clear figure(1)
save(color_vec_path, 'color_vec', '-ASCII')
close all