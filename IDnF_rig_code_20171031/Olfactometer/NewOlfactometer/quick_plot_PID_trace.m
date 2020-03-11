path = 'E:\Turner lab\Data\mehrab\20200203\handover_PID_traces_set2\';

a = get_PID_traces(path, 55, 0.099);
a(:, 1) = a(:, 1) - a(:, 2);

plot(a, 'lineWidth', 2)