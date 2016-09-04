clear all
close all

frac_vec = 0.0:0.1:1
for active_frac_iter = 1:10
    active_frac = frac_vec(active_frac_iter);
    n_iters = 100;
    saved_data = zeros(3, n_iters);
    for r_iter = 1:n_iters

        set1 = randperm(100);
        set1 = set1(1:round(active_frac.*100));

        set2 = randperm(100);
        set2 = set2(1:round(active_frac.*100));

        overlap = length(intersect(set1, set2));
        nonoverlap = (length(set1).*2) - overlap;

        ratio = overlap./nonoverlap;

        saved_data(:, r_iter) = [overlap; nonoverlap; ratio];

    end

    saved_data_main(:, active_frac_iter) = [mean(saved_data(1, :) ), mean(saved_data(2, :) ), mean(saved_data(3, :) )];
end
plot(saved_data_main(3, :))