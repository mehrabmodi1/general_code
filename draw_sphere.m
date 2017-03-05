function [] = draw_sphere(fig_n, center_trip, radius, color_val, alpha_val)

[x, y, z] = sphere;
c_map = zeros(size(z, 1), size(z, 2));
c_map(:,:) = color_val;
c_map(1,1) = 0;
c_map(1,2) = 1;

figure(fig_n)
surf((radius.*x + center_trip(1)), (radius.*y + center_trip(2)), (radius.*z + center_trip(3)), c_map, 'edgeColor', 'none');
alpha(alpha_val)