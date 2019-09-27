direc = 'C:\Users\Mehrab\Dropbox (HHMI)\data_sharing\Data_for_Herve\KC_AlphaBeta_PABAEL_201908set\fly3\';

data = load([direc, 'dFF_data.mat']);

dff_data_mat_f = data.dff_data_mat_f;
stim_mat = load([direc, 'stim_mat.mat']);
stim_mat = stim_mat.stim_mat;

a = dff_data_mat_f(:, :, 1);
b = dff_data_mat_f(:, :, 5);

difsum = sum(sum(a - b, 'omitnan'),'omitnan')