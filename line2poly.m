function [poly_crds] = line2poly(x_crds, y_crds, thickness)
%This function takes a line specified in 2D by its x-coordinates and
%y-coordinates and generates a polygon that traces the path of the line at
%a fixed, specified thickness.
%Mehrab N. Modi 20131011


% x_crds = [1:1:10];
% y_crds = [2, 4, 6, 3, 1, 5, 9, 9, 10, 0];
% thickness = 0.1;
% 
% 
% figure(1)
% plot(x_crds, y_crds)

% 
% if nargin < 2
%     error('line must be specified in 2 D')
% elseif nargin < 3
%     disp('thickness assumed to be 1')
%     thickness = 1;
% else
% end

small_step_x = mean(diff(x_crds))./500;     %calculates a negligibly small step-size along x


%correcting for non 1:1 scale on x and y by normalising x_crds and y_crds now, and reversing this at the end
x_max = nanmax(x_crds);
y_max = nanmax(y_crds);
x_crds = x_crds./x_max;
y_crds = y_crds./y_max;


polygon_points1 = [];
polygon_points2 = [];
for point_no = 1:(length(x_crds)-1)
    %calculating angle between line joining curr_point and next point and
    %horizontal
    curr_point = [x_crds(point_no), y_crds(point_no)];
    next_point = [x_crds(point_no + 1), y_crds(point_no + 1)];
    
    if max(isnan(curr_point)) == 1
        continue
    elseif max(isnan(curr_point)) == 1
        continue
    else
    end
    
    curr_vec = next_point - curr_point;     %subtracting away co-ords of first point to make it the origin
    ref_vec = [1, 0];                       % reference, horizontal vector
    Cos_h_angle = dot(curr_vec, ref_vec)/(norm(curr_vec)*norm(ref_vec));
    h_angle = acos(Cos_h_angle).*(180/pi);             %angle of current segment rel to the horizontal in degrees
      

    %for segment to have desired thickness, distance along y between current pair of polygon vertices must be adjusted.
    y_sep = thickness./sind(90 - h_angle);      %adjusted y-separation to maintain required line thickness
    
    %building polygon - assigning 4 vertices for each point specified by
    %x_crds and y_crds
    
    %adding co-ords for first point in current line-segmen  t
    polygon_points1 = [polygon_points1; (curr_point(1, 1) + small_step_x), curr_point(1, 2) + (y_sep./2)];      %set of points above specified line
    polygon_points2 = [polygon_points2; (curr_point(1, 1) + small_step_x), curr_point(1, 2) - (y_sep./2)];      %set of points below specified line
    %adding co-ords for second point in current line-segment
    polygon_points1 = [polygon_points1; (next_point(1, 1) - small_step_x), next_point(1, 2) + (y_sep./2)];      %set of points above specified line
    polygon_points2 = [polygon_points2; (next_point(1, 1) - small_step_x), next_point(1, 2) - (y_sep./2)];      %set of points below specified line
    
    
end

poly_crds = [polygon_points1; flipud(polygon_points2)];

%reversing initial normalisation
poly_crds(:, 1) = poly_crds(:, 1).*x_max;
poly_crds(:, 2) = poly_crds(:, 2).*y_max;

% figure(2)
% patch(poly_crds(:, 1), poly_crds(:, 2), 'r')
%         