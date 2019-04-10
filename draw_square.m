function [im] = draw_square(im, x, y, size, fill)
%This script draws a square of side length = 2*size + 1, centered at x, y
%in the matrix im, replacing these values with the value fill. 
%Mehrab Modi, 20141104

x = round(x);
y = round(y);

for x_pix = (-1.*size):size
    for y_pix = (-1.*size):size
        im( (y + y_pix), (x + x_pix) ) = fill;        
    end
end


end