clear all
close all

path = 'C:\Data\Data\Analysed_data\Suite2P_results\20180514\odpulse_freq_series\1\odor_trs_00003.tif';

frame_orig = imread(path, 1);

frame = im2double(frame_orig);