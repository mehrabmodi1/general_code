function [] = plot_session_resp_trend(resp_mat, paired_tr_n, color_vecs, fig_n)
%syntax: [] = plot_session_resp_trend(resp_mat, paired_tr_n, color_vecs, fig_n)
%This function plots response sizes averaged across flies, trial by trial
%as imaging sessions proceed.

%resp_mat description
%dim1 - repeats
%dim2 - odor1pre, odor1post, odor2pre...
%dim3 - fly_n

%averaging each repeat across flies
mean_mat = mean(resp_mat, 3, 'omitnan');
se_mat = std(resp_mat, [], 3, 'omitnan')./sqrt(size(resp_mat, 3));

%concatenating column pairs for post resps on at the end of the pre resps
pre_col_vec = 1:2:size(mean_mat, 2);
post_col_vec = pre_col_vec + 1;
mean_mat_reshaped = cat(1, mean_mat(:, pre_col_vec), mean_mat(:, post_col_vec));
se_mat_reshaped = cat(1, se_mat(:, pre_col_vec), se_mat(:, post_col_vec));

%identifying and removing pad for any trials not acquired for all odors
temp = mean(mean_mat_reshaped, 2, 'omitnan');
rem_trs = find(isnan(temp) == 1);
mean_mat_reshaped(rem_trs, :) = [];
se_mat_reshaped(rem_trs, :) = [];

%plotting response trends
figure(fig_n)
for od_n = 1:size(mean_mat_reshaped, 2)
    shadedErrorBar([], mean_mat_reshaped(:, od_n), se_mat_reshaped(:, od_n), {'Color', color_vecs(od_n, :)}, 1)
    hold on
end

ylabel('mean dF/F response')
xlabel('repeat number')
ax_vals = axis;
ax_vals(1, 3) = min([0, ax_vals(1, 3)]);
axis(ax_vals)
%plotting pairing trial indicator
plot([paired_tr_n, paired_tr_n], [ax_vals(3), ax_vals(4)], '--r', 'lineWidth', 2)
