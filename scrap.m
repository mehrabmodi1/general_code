clear all
close all

%color_vec = load('C:\Data\Code\general_code\std_color_vec.txt');
color_vec = [0.5529    0.8275    0.7804;
                     0.8500    0.8500    0.6520;
                     0.4    0.4    0.4;
                     0.9843    0.5020    0.4471;
                     0.5020    0.6941    0.8275;
                     0.9922    0.7059    0.3843;
                     0.7020    0.8706    0.4118;
                     0.9882    0.8039    0.8980;
                     0.8510    0.8510    0.8510;
                     0.7373    0.5020    0.7412;
                     0.8000    0.9216    0.7725;
                     0.8750    0.8294    0.4353;];
%color_vec = color_vec./255;
n_colors = 12;
rand_mat = rand(100, n_colors);
add_mat = repmat( (1:n_colors), 100, 1);
rand_mat = rand_mat + add_mat;

for line_n = 1:n_colors
    plot(rand_mat(:, line_n), 'Color', color_vec(line_n, :), 'LineWidth', 3)
    hold on
end
hold off