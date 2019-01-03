clear all
close all
figure(1)
[del, img1] = vennX_MM( [ 20 10 20 10 20 10 10], .01, [], 0, 1, 1);
alpha(0.5)
hold on

[del, img2] = vennX_MM( [ 25 10 15 10 20 10 10], .01, [], 0, 1, 1 );
alpha(0.5)

figure(2)
a = imagesc(img1);
hold on
b = imagesc(img2);
alpha(b, 0.5)
