function [z_stack_mat] = load_z_stack()

    [del, m_aq_direc] = curr_aq_direc;
    if exist([m_aq_direc, 'zstack\z_stack.mat']) == 0
        trial_n = 1;
        dir_contents = dir([m_aq_direc, 'zstack\*.tif']);
        tif_info = imfinfo([m_aq_direc, 'zstack\', dir_contents(1).name]);
        keyboard
        n_planes = size(tif_info, 1);

        plane = imread([m_aq_direc, 'zstack\', dir_contents(1).name], 1);
        z_stack_mat = zeros(size(plane, 1), size(plane, 2), n_planes);
        z_stack_mat(:, :, 1) = plane;

        for plane_n = 2:n_planes
            plane = imread([m_aq_direc, 'zstack\', dir_contents(1).name], plane_n);
            z_stack_mat(:, :, plane_n) = plane;
        end
        save([m_aq_direc, 'zstack\z_stack.mat'], 'z_stack_mat');
    
    elseif exist([m_aq_direc, 'zstack\z_stack.mat']) == 2
        trial_n = 2;
        z_stack_mat = load([m_aq_direc, 'zstack\z_stack.mat']);
        z_stack_mat = z_stack_mat.z_stack_mat;
    end
    clear plane
    clear del



end