function [score_vec] = compare_with_zstack(curr_plane, z_stack, plot_scores, alg_n)
%This function compares a given slice with the entire z-stack and returns
%the best matching slice number from the z stack.

    n_planes = size(z_stack, 3);
    score_vec = zeros(n_planes, 1);
    curr_plane(curr_plane <= 2) = 0;
    for plane_n = 1:n_planes
        norm_plane = z_stack(:, :, plane_n);
        norm_plane(norm_plane <= 2) = 0;
        
        figure(999)
        subplot(2, 1, 1)
        imagesc(curr_plane)
        subplot(2, 1, 2)
        imagesc(norm_plane)
        
        del = input('press enter');
        
        if alg_n == 1
            corr_mat = xcorr2_fft(curr_plane, norm_plane);
            score_vec(plane_n, 1) = max(max(corr_mat));
        elseif alg_n == 2
            corr_val = mat_corrcoef(curr_plane, norm_plane);
            score_vec(plane_n, 1) = corr_val;
        else
        end
    end
if plot_scores == 1
    plot(score_vec)
else
end

end