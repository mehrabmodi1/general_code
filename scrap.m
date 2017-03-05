clear all
close all

black_dist = normrnd(5, 2, 5000, 1);
grey_dist = normrnd(1, 2, 5000, 1);

figure(1)
histogram(black_dist)
hold on
histogram(grey_dist)

thresh1 = 3.5;
thresh2 = .2;

%using higher threshold first (same as in paper)
black_high1 = length(find((black_dist + .2) > thresh1));
grey_high1 = length(find((grey_dist + .2) > thresh1));

black_high2 = length(find(black_dist > thresh1));
grey_high2 = length(find(grey_dist > thresh1));

disp('difference between counts above higher thresh; with offset added')
diff1 = black_high1 - grey_high1

disp('difference between counts above higher thresh; without offset added')
diff2 = black_high2 - grey_high2



%---------------
%lower threshold (lower than mean of grey dist)
black_low1 = length(find((black_dist + .2) > thresh2));
grey_low1 = length(find((grey_dist + .2) > thresh2));

black_low2 = length(find(black_dist > thresh2));
grey_low2 = length(find(grey_dist > thresh2));

disp('difference between counts above lower thresh; with offset added')
diff1 = black_low1 - grey_low1
disp('difference between counts above lower thresh; without offset added')
diff2 = black_low2 - grey_low2

