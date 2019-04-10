function [] = cloud_plotter(fig_n, data_mat, cloud_dias, alpha_val)

n_clouds = size(data_mat, 2);
color_vec = 0:(1./n_clouds):1;

for cloud_n = 1:n_clouds
    curr_color = color_vec(cloud_n);
    curr_dia = cloud_dias(cloud_n);
    figure(fig_n)
    draw_sphere(fig_n, [data_mat(:, cloud_n)], curr_dia, curr_color, alpha_val);
    hold on
    %keyboard
end
grid off
hold off