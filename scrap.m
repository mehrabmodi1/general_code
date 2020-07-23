im = rand(100, 100);
figure(3)
bw = roipoly(im);
figure(4)
x = [10, 50, 50, 10];
y = [10, 10, 50, 50];

bw2 = roipoly(im, x, y);
%bw3 = poly2mask(x, y, 100, 100);
figure(5)
imagesc(bw2)