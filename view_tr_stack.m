function [] = view_tr_stack(direc)

tr_stack = load([direc, 'tr_avg_stack.mat']);
tr_stack = tr_stack.ave_stack;

playStack(tr_stack, 0.2, []);
