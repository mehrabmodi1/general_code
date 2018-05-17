clear all
close all

path = 'E:\Data\Raw_Data_Current\Resonant\20180515\OK107_x_jGCaMP7f_freq_series\1\odor_trs_00001.tif';

stack_obj = ScanImageTiffReader('E:\Data\Raw_Data_Current\Resonant\20180515\OK107_x_jGCaMP7f_freq_series\1\odor_trs_00001.tif');
stack = stack_obj.data();
stack = permute(stack,[2 1 3]);

path2 = 'E:\Data\Analysed_data\Suite2p\Results\20180515\OK107_x_jGCaMP7f_freq_series\1\F_20180515_OK107_x_jGCaMP7f_freq_series_plane1_proc';
