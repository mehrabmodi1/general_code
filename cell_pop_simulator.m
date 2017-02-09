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

%%
%RE-SAMPLING
%generating a re-sampled population response dataset
direc_n = 1;
n_cells = 300;
n_odors = 100;
m_sparseness = .1;       %range: 0 to 1          defined as a property of each odor
sd_sparseness = .05;      %range: 0 to 1          defined as a property of each odor

%This parameter controls how likely a particular cell is to respond to more than one odor
cooperativity = -1;        %range: -1 to 1        %for  0 to 1; defines the fraction of sig odor resps that are picked up from perfect cooperativity and re-distributed randomly
                                                %for -1 to 0; defines the fraction of sig odor resps that are picked up from a systematically distributed set of responses and re-distributed randomly
                                                
%function to generate re-sampled population responses
sim_data_mat = resample_cell_pop(orig_data_mat{direc_n, 1}, n_cells, n_odors, m_sparseness, sd_sparseness, cooperativity);
keyboard
%function to take a re-sampled population response matrix and do a PCA and calculate odor separablity metrics
%classification_scores = odor_classifier(sim_data_mat);



