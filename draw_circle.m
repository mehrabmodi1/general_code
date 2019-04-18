function [matrix_out, points_list] = draw_circle(centrex, centrey, radius, matrix, filled)
%Mehrab Modi - 3/9/13 syntax: [matrix_out, points_list] = draw_circle(centrex, centrey, radius, matrix, filled)
%this function draws a circle in a given matrix, with a given centre and
%radius. it also gives as a second output, a list of the points that fall
%within the circle, or on its perimeter - this last factor being determined
%by the 'filled' input, with 1 meaning a filled circle and 0 just an
%outline.

if nargin < 4
    disp('not enough inputs')
    
elseif nargin == 4
    disp('circle will be filled')
    filled = 1;
else
end

sizex = size(matrix, 1);
sizey = size(matrix, 2);

centrex = floor(centrex);
centrey = floor(centrey);

%building list of possible x co-ordinate values
if (centrex + radius) < sizex && (centrex - radius) > 1
    x_min = centrex - radius;
    x_max = centrex + radius;
elseif (centrex - radius) < 2 && (centrex + radius) < sizex
    disp('circle crossing matrix edge, clipped circle drawn.')
    x_min = 1;
    x_max = centrex + radius;
elseif (centrex + radius) > sizex && (centrex - radius) > 1
    disp('circle crossing matrix edge, clipped circle drawn.')
    x_min = centrex - radius;
    x_max = sizex;
elseif (centrex + radius) > sizex && (centrex - radius) < 2
    x_min = 1;
    x_max = sizex;
else
end



 %building list of possible y co-ordinate values
if (centrey + radius) < sizey && (centrey - radius) > 1
    y_min = centrey - radius;
    y_max = centrey + radius;
elseif (centrey - radius) < 2 && (centrey + radius) < sizey
    disp('circle crossing matrix edge, clipped circle drawn.')
    y_min = 1;
    y_max = centrey + radius;
elseif (centrey + radius) > sizey && (centrey - radius) > 1
    disp('circle crossing matrix edge, clipped circle drawn.')
    y_min = centrey - radius;
    y_max = sizey;
elseif (centrey + radius) > sizey && (centrey - radius) < 2
    y_min = 1;
    y_max = sizey;
else
end  
    

%using equation of circle to find points on circle
points_list = [];
c_point_no = 0;
for x = x_min:1:x_max
    for y = y_min:1:y_max
        if filled == 1
            if ( (x - centrex)^2 + (y - centrey)^2) <= radius^2;
                c_point_no = c_point_no + 1;
                points_list = [points_list; x, y]; 
            else
            end
        elseif filled == 0
            if round((x - centrex)^2 + (y - centrey)^2) == radius^2;
                c_point_no = c_point_no + 1;
                points_list = [points_list; x, y]; 
            else
            end
        else
        end
    end
end




matrix_out = matrix;

for pt_no = 1:length(points_list)
    
    pt = points_list(pt_no, :);
    
    matrix_out(pt(1, 1), pt(1, 2)) = 1;
end