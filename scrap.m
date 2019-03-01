clear all
close all

f = figure;
axis([0, 10, 0, 5])

greenrect = rectangle('Position',[0.7 0.1 0.1 0.2]);
set(greenrect,'FaceColor','g');
draggable(greenrect);

hold on

a = rand(2, 2) + 1;
im = imagesc(a);
draggable(im, 'none', [-inf inf -inf inf]);

im.XData
im.YData