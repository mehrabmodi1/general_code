clear all
close all

%model parameters
%n_odors
%n_cells

%sparseness_mn   range: 0 to 1          defined as a property of each odor
%sparseness_sd   range: 0 to 1          defined as a property of each odor

%odor_resp_overlap

%list of data dump direcs to load cell data from
dump_direcs = {'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Analysis\data_dump\AB';... %AB cells
               'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Analysis\data_dump\ApBp';... %ApBp cells
               'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Analysis\data_dump\G';... %G cells
               'C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Analysis\data_dump\OK107'... %OK107 cells
                };

%loading orig cell data
for dump_direc_n = 1:length(dump_direcs)
    curr_dir = dump_direcs{dump_direc_n, 1};
    data_mat = load([curr_dir '\cell_data.mat']);
    data_mat = data_mat.saved_data;
    orig_data_mat{dump_direc_n, 1} = data_mat;                            %mat in which all orig data will be stored in memory
    
end
    
%function to generate re-sampled population responses
%PICK UP THREAD HERE
%write function to generate resampled pop resps
