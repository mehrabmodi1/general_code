clear all
close all
path1 = 'C:\Users\Mehrab\Dropbox (HHMI)\data_sharing\Data_for_Herve\discrimination_data_Herve\Gamma\fly3\';
path2 = 'C:\Users\Mehrab\Dropbox (HHMI)\data_sharing\Data_for_Herve\KC_Gamma_PABAEL_201908set\fly1\';

stim_mat1 = load([path1, 'stim_mat.mat']);
stim_mat1 = stim_mat1.stim_mat;
stim_mat2 = load([path2, 'stim_mat.mat']);
stim_mat2 = stim_mat2.stim_mat;
stim_mat_simple_2 = generate_stim_mat_simple(stim_mat2);
