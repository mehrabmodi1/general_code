function [stack_mat] = load_tif_stack(path)
    tif_info = imfinfo(path);
    n_planes = size(tif_info, 1);

    plane = imread(path, 1);
    stack_mat = zeros(size(plane, 1), size(plane, 2), n_planes);
    stack_mat(:, :, 1) = plane;

    for plane_n = 2:n_planes
        plane = imread(path, plane_n);
        stack_mat(:, :, plane_n) = plane;
    end


end
