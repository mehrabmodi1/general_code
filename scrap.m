clear all
close all

x = zeros(20, 20);
x1 = [3, 21];
x2 = [15, 15];
interp_vals = round([linspace(x1(1), x2(1), 20); linspace(x1(2), x2(2), 20)]);

x(x1(1), x1(2)) = 1;
x(x2(1), x2(2)) = 1;

y = zeros(20, 20);

for i_val = 1:size(interp_vals, 2)
    y(interp_vals(1, i_val), interp_vals(2, i_val)) = 1;
    
end

figure(1)
imagesc(x)

figure(2)
imagesc(y)