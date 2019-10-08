function [ROI_groups] = define_roi_groups(ROI_mat, ave_im, ignore_ROIs)
%This function pops up a figure showing all the ROIs in ROI_mat, and asks
%the user to click on and define ROIs that belong in a group together. The
%user decides how many ROIs go in each group, and how many groups there
%are.
% clear all
% close all
% 
% ave_im = rand(1000, 1000);
% ROI_mat = zeros(1000, 1000, 5);
% ROI_mat(200:300, 200:300, 1) = 1;
% ROI_mat(200:300, 700:900, 2) = 1;
% ROI_mat(500:600, 200:300, 3) = 1;
% ROI_mat(200:300, 500:600, 4) = 1;
% ROI_mat(700:900, 200:300, 5) = 1;

%creating manual selection figure
ROI_mat_summed = sum(ROI_mat, 3);
cyan_im = cat(3, zeros(size(ave_im)), ones(size(ave_im)), (zeros(size(ave_im)) + 0.5));
figure(1)
imshow(imadjust(ave_im))
hold on
h = imshow(cyan_im);
hold off
set(h, 'AlphaData', ROI_mat_summed.*0.4) 
%set(1, 'Position', [200, 50, 100 + 800, 100 + 800]);
plot_big_fig(1)

group_n = 1;
prev_gp_n = 1;
assigned_ROIs = [];
counter = 0;
while 1
    counter = counter + 1;
    disp(['current group number: ' int2str(group_n)]);
    curr_pos = ginput(1);
    curr_pos = round(curr_pos);
    
    if group_n ~= prev_gp_n
        ROI_groups{group_n, 1} = [];
    elseif counter == 1
        ROI_groups{group_n} = [];
    else
    end
    prev_gp_n = group_n;
    
    ROI_set = 1:1:size(ROI_mat, 3);
    ROI_set(ignore_ROIs) = [];
   for ROI_ni = 1:size(ROI_set, 1)
       ROI_n = ROI_set(ROI_ni);
        if ROI_mat(curr_pos(2), curr_pos(1), ROI_n) == 1
            %removing ROI from any groups it has already been assigned to
            for group_n_temp = 1:size(ROI_groups, 1)
                curr_mems_temp = ROI_groups{group_n_temp};
                del = find(curr_mems_temp == ROI_n);
                curr_mems_temp(del) = [];
                ROI_groups{group_n_temp, 1} = curr_mems_temp;
            end
            
            %assigning currently selected ROI to current group
            ROI_groups{group_n, 1} = [ROI_groups{group_n}, ROI_n];
            if isempty(intersect(assigned_ROIs, ROI_n)) == 1
                assigned_ROIs = [assigned_ROIs; ROI_n];
            else
            end
        else
        end
    end
    
    assigned_ROIs_temp = assigned_ROIs;
    curr_mems = ROI_groups{group_n, 1};
    [del, to_remove] = intersect(assigned_ROIs_temp, curr_mems);
    assigned_ROIs_temp(to_remove) = [];                     %this vector lists ROIs that have been assigned, but not to the current group
    
    other_ROIs = 1:1:size(ROI_mat, 3);
    other_ROIs(assigned_ROIs) = [];                         %this vector lists un-assigned ROIs
    
    ROI_mat_summed_bk = sum(ROI_mat(:, :, other_ROIs), 3);
    ROI_mat_summed_assigned = sum(ROI_mat(:, :, assigned_ROIs_temp), 3);
    ROI_mat_summed_curr = sum(ROI_mat(:, :, curr_mems), 3);
    
    cyan_im = cat(3, zeros(size(ave_im)), ones(size(ave_im)), (zeros(size(ave_im)) + 0.5));
    figure(1)
    imshow(imadjust(ave_im))
    hold on
    h = imshow(cyan_im);
    hold off
    set(h, 'AlphaData', ROI_mat_summed_bk.*0.4) 
    %set(1, 'Position', [200, 50, 100 + 800, 100 + 800]);
    plot_big_fig(1)

    hold on
    red_im = cat(3, (zeros(size(ave_im)) + 1), zeros(size(ave_im)), (zeros(size(ave_im)) + 0.5));
    h = imshow(red_im);
    hold off
    set(h, 'AlphaData', ROI_mat_summed_curr.*0.4) 
    %set(1, 'Position', [200, 50, 100 + 800, 100 + 800]);
    plot_big_fig(1)

    hold on
    yellow_im = cat(3, (zeros(size(ave_im)) + 1), (zeros(size(ave_im)) + 0.7), (zeros(size(ave_im)) ));
    h = imshow(yellow_im);
    hold off
    set(h, 'AlphaData', ROI_mat_summed_assigned.*0.4) 
    %set(1, 'Position', [200, 50, 100 + 800, 100 + 800]);
    plot_big_fig(1)
    
    %ending loop if all ROIs have been assigned to groups.
    if length(assigned_ROIs) == size(ROI_mat, 3)
        break
    else
    end
    
    button = questdlg('What would you like to do?','Grouping ROIs','Continue grp','New grp','End grouping', 'Continue grp');
    
    if strcmp(button, 'End grouping') == 1
        break
    elseif strcmp(button, 'Continue grp') == 1
        group_n = group_n;
    elseif strcmp(button, 'New grp') == 1
        group_n = group_n + 1;
    else
    end
    
end
close figure 1 
%end