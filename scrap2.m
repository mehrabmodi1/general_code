clear all
close all

base_path1 = 'C:\Data\Data\Analysed_data\data_sharing\Rob_data_fromShyam\campbell_data_20220627\';
file1 = 'TSeries-01052011-1157-001.mat';

base_path2 = 'C:\Data\Data\Analysed_data\data_sharing\Rob_data_fromShyam\campbell_data_20220627\TSeries-01052011-1157-001\';
file2 = 'params_110107_152032.mat';

data_struc = load([base_path2, file2]);
%data_struc = data_struc.data;