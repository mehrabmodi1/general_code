clear all
close all

global xpos1
global ypos1
global xpos0
global ypos0

f = figure;
bk_mat = rand(100, 100);
imagesc(bk_mat);
colormap('gray')

hold on

a = zeros(10, 10) + 2;
im = imagesc(a);
alpha(im, 0.5)
draggable(im, 'none', [-inf inf -inf inf], 'endfcn', @get_final_pos);

xpos0 = im.XData;
ypos0 = im.YData;

function get_final_pos(im)
    disp('endfcn executed')
    global xpos1
    xpos1 = im.XData;
    global ypos1 
    ypos1 = im.YData;
    global xpos0
    global ypos0
    
    ydisp = ypos1 - ypos0
    xdisp = xpos1 - xpos0
end