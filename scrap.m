clear all
close all
<<<<<<< Updated upstream
path1 = 'C:\Users\Mehrab\Dropbox (HHMI)\data_sharing\Data_for_Herve\discrimination_data_Herve\Gamma\fly3\';
path2 = 'C:\Users\Mehrab\Dropbox (HHMI)\data_sharing\Data_for_Herve\KC_Gamma_PABAEL_201908set\fly1\';

stim_mat1 = load([path1, 'stim_mat.mat']);
stim_mat1 = stim_mat1.stim_mat;
stim_mat2 = load([path2, 'stim_mat.mat']);
stim_mat2 = stim_mat2.stim_mat;
stim_mat_simple_2 = generate_stim_mat_simple(stim_mat2);
=======

path = 'E:\Data\Raw_Data_Current\Resonant\20200821\fly1_d5HT1BKC_handover_rtrains_prehabit\1\';

stim_mat = load([path, 'params.mat']);
stim_mat = stim_mat.params_mat;

stim_mat_simple = generate_stim_mat_simple(stim_mat);
>>>>>>> Stashed changes
