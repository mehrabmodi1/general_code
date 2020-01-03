function [im_warped, im_orig, landmark_pairs] = warp_im(im_orig, landmark_pairs)
%This function takes an image and a set of pairs of points (landmark_pairs, size n, 4)
%identified in the original image and a second, warped reference frame of
%the same size as the original image. It then warps each point in the
%original image to bring each landmark point close to it's new location in
%the warped reference frame. Each pixel is shifted by the mean
%displacements of all the landmarks, weighted by distance from the pixel
%being shifted. Unassigned pixels in the warped image are interpolated from
%neighbouring pixels.

%Test inputs
% im_orig = zeros(100, 100) + nan;
% im_orig(10:20, 10:20) = 1;
% im_orig(10:20, 80:90) = 1;
% im_orig(80:90, 10:20) = 1;
% im_orig(80:90, 80:90) = 1;
% landmark_pairs(1, :) = [5, 5, 20, 20];
% landmark_pairs(2, :) = [85, 85, 50, 50];


%this sets the exp decay rate of the weighting of landmarks by distance
one_eth_decay_dist = 0.2;

%computing the displacement of each landmark point
lm_displacements = [(landmark_pairs(:, 1) - landmark_pairs(:, 3)), (landmark_pairs(:, 2) - landmark_pairs(:, 4))];
max_dist = round(sqrt(size(im_orig, 1).^2 + size(im_orig, 2).^2));  %the length of the diagonal of im orig
decay_rate = 1./(max_dist.*one_eth_decay_dist);

%computing new locations for each pixel and assigning pixel values in
%warped image
im_warped = zeros(size(im_orig, 1), size(im_orig, 2)) + nan;
im_warped_pix_counts = zeros(size(im_orig, 1), size(im_orig, 2));       %to keep track of how many pixels' values from im_orig have been assigned to each pixel in im_warped
for row_n = 1:size(im_orig, 1)
    for col_n = 1:size(im_orig, 2)
        coords = [row_n, col_n; landmark_pairs(:, 1:2)];
        dist_vec = squareform(pdist(coords));
        dist_vec = dist_vec(1, :);
        dist_vec(1) = [];               %vector of distances of curr point in im_orig from each landmark in im_orig
        landmark_wts = 1./exp(decay_rate.*dist_vec);
        
        %weighted mean displacement for current pixel in im_orig to be inserted into im_warped
        curr_displacements = round(weighted_mean(lm_displacements, landmark_wts'));
        
        %inserting current pixel from im_orig into im_warped if possible.
        %If another pixel has already been inserted at the same location in
        %im_warped, the mean of the two will be taken.
        new_pos = [row_n - curr_displacements(1), col_n - curr_displacements(2)];
        
        %checking if new position falls outside im_warped and skipping
        if min(sign(new_pos)) <= 0
            continue
        elseif new_pos(1) > size(im_orig, 1) || new_pos(2) > size(im_orig, 2)
            continue
        else
        end
        
        %assigning value
        im_warped(new_pos(1), new_pos(2)) = sum([im_orig(row_n, col_n), im_warped(new_pos(1), new_pos(2))], [], 'omitnan');
        im_warped_pix_counts(new_pos(1), new_pos(2)) = im_warped_pix_counts(new_pos(1), new_pos(2)) + 1;
        
%         if row_n == 85 && col_n == 85
%             keyboard
%         else
%         end
%         
    end
end
im_warped = im_warped./im_warped_pix_counts;