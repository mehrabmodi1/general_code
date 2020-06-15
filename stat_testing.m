
function [p_vec_out, effect_sizes_out, min_dif_out] = stat_testing(data_mat, is_paired, non_param, p_val)
%Syntax: [p_vec_out, effect_sizes_out, min_dif_out] = stat_testing(data_mat, is_paired, non_param, p_val)
%This function expects data in pairs of columns to be compared. It tests
%for paired or independent samples based on the boolean is_paired, and uses
%normal or non-parametric tests based on the boolean non_param. It computes
%needed n based on the maximum p_val input, for norm tests. p_vec_out and effect_sizes out are
%vectors of p values and Cohen's d of length n column pairs. min_dif_out is
%a matrix of size n column pairs x 4 ie. [min difference, actual
%difference, needed n, actual n]. For non parametric tests, only p_vec_out
%is defined.

    n_col_pairs = size(data_mat, 2)./2;
    p_vec_out = zeros(n_col_pairs, 1) + nan;
    effect_sizes_out = zeros(n_col_pairs, 1) + nan;
    for col_pair_n = 0:(n_col_pairs - 1)
        col1 = (col_pair_n.*2) + 1;
        col2 = col1 + 1;
            n = size(data_mat, 1);
            
            if non_param == 0   %means, assume normal distributions
                if is_paired == 1       %paired samples
                    [del, p] = ttest(data_mat(:, col1), data_mat(:, col2));
                    effect_size = computeCohen_d(data_mat(:, col1), data_mat(:, col2), 'paired');
                    paired_difs = data_mat(:, col1) - data_mat(:, col2);    %for paired-samples, sample distribution is distribution of paired differences.
                    mean1 = mean( paired_difs, 'omitnan');                  %actual mean of paired differences
                    sd1 = std(paired_difs, [], 'omitnan');
                    min_mean1 = sampsizepwr('t',[0 sd1],[],(1 - p_val.*2), n);        %null hypothesis mean is difference of 0
                    needed_n = sampsizepwr('t',[0 sd1], mean1, (1 - p_val.*2), []);
                    min_dif = abs(min_mean1 - 0); 
                    ac_dif = mean1;
                elseif is_paired == 0   %independent samples
                    [del, p] = ttest2(data_mat(:, col1), data_mat(:, col2));
                    effect_size = computeCohen_d(data_mat(:, col1), data_mat(:, col2), 'independent');
                    mean1 = mean(data_mat(:, col1), 'omitnan');
                    ac_mean2 = mean(data_mat(:, col2), 'omitnan');
                    sd1 = std(data_mat(:, col1), [], 'omitnan');
                    mean2 = sampsizepwr('t2',[mean1 sd1],[],(1 - p_val.*2), n);        %two-tailed test
                    needed_n = sampsizepwr('t',[0 sd1], mean1, (1 - p_val.*2), []);
                    min_dif = abs(mean2 - mean1);
                    ac_dif = abs(ac_mean2 - mean1);
                else
                end
            elseif non_param == 1   %medians, non-parametric tests
                if is_paired == 1       %paired samples
                    p = signrank(data_mat(:, col1), data_mat(:, col2));
                    effect_size = [];
                    min_dif = [];
                    ac_dif = [];
                    needed_n = [];
                    n = [];
                elseif is_paired == 0   %independent samples
                    p = ranksum(data_mat(:, col1), data_mat(:, col2));    %Mann Whitney U test
                    effect_size = [];
                    min_dif = [];
                    needed_n = [];
                    ac_dif = [];
                    n = [];
                else
                end
                
            else
            end
        try
            p_vec_out((col_pair_n + 1), 1) = p;
            effect_sizes_out((col_pair_n + 1), 1) = effect_size;
            min_dif_out((col_pair_n + 1), 1:4) = [min_dif, ac_dif, needed_n, n];
        catch
            keyboard
        end
    

    end
    
    if is_paired == 1   %paired-samples
        if non_param == 0   %use means, parametric test
            disp('Mean with SEs, paired-sample T test.')
        elseif non_param == 1  %use medians, non-parametric test
            disp('Median with quartiles, Wilcoxon signed rank test.')
        else
        end
    elseif is_paired == 0   %independent samples
        if non_param == 0   %use means, parametric test
            disp('Mean with SEs, two-sample T test.')
        elseif non_param == 1  %use medians, non-parametric test
            disp('Median with quartiles, Mann Whitney U test.')
        else
        end
    else
    end
        

end