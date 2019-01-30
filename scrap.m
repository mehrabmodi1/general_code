naive = [0.6, 1, 0.8, 0.7, 0.9, 0.7, 0.9, 0.9, 1, 0.8];
naive = naive';
vecs = repmat(naive, 1, 8);

vecs(2, 2) = 0.7;
vecs(7, 2) = 0.6;
vecs(10, 2) = 0.6;

vecs(2, 3) = 0.7;
vecs(7, 3) = 0.6;
vecs(10, 3) = 0.6;
vecs(3, 3) = 0.9;
vecs(5, 3) = 1;
vecs(9, 3) = 1;

vecs(2, 4) = 0.5;
vecs(7, 4) = 0.4;
vecs(10, 4) = 0.4;
vecs(3, 4) = 0.9;
vecs(5, 4) = 1;
vecs(9, 4) = 1;

vecs(3, 5) = 1;
vecs(5, 5) = 1;
vecs(9, 5) = 1;
vecs(2, 5) = 0.5;
vecs(7, 5) = 0.4;
vecs(10, 5) = 0.4;

vecs(2, 6) = 0.4;
vecs(7, 6) = 0.3;
vecs(10, 6) = 0.3;
vecs(3, 6) = 1;
vecs(5, 6) = 1;
vecs(9, 6) = 1;

vecs(2, 7) = 0.4;
vecs(7, 7) = 0.3;
vecs(10, 7) = 0.3;
vecs(3, 7) = 1;
vecs(5, 7) = 1;
vecs(9, 7) = 1;

vecs(2, 8) = 0.0;
vecs(7, 8) = 0.0;
vecs(10, 8) = 0.0;
vecs(3, 8) = 1;
vecs(5, 8) = 1;
vecs(9, 8) = 1;

a = pdist(vecs');
a = squareform(a);

figure(1)
plot(a(:, 1), 'k', 'lineWidth', 3)
fig_wrapup(1)
