mat_ser = load('E:\Data\Analysed_data\Suite2p\Results\20180515\OK107_x_jGCaMP7f_freq_series\1\extracted_raw_data_mat_old');
mat_ser = mat_ser.raw_data_mat;
mat_par = load('E:\Data\Analysed_data\Suite2p\Results\20180515\OK107_x_jGCaMP7f_freq_series\1\extracted_raw_data_mat');
mat_par = mat_par.raw_data_mat;
 
diff_mat = mat_ser(1:size(mat_par, 1), :, 1) - mat_par(:, :, 1);

tot_diff = sum(sum(diff_mat))

trace_ser = mat_ser(1:size(mat_par, 1), 1, 1);
trace_par = mat_par(:, 1, 1);

figure(1)
plot(trace_ser)
hold on
plot(trace_par + std(trace_ser))

