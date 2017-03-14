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
n_cells = 100;
n_odors = 50;
m_sparseness = .2;        %range: 0 to 1          sparseness of each odor drawn from this dist
sd_sparseness = .05;      %range: 0 to 1          sparseness of each odor drawn from this dist

%This parameter controls how likely a particular cell is to respond to more than one odor
cooperativity = 0;       %range: -1 to 1        %for  0 to 1; defines the fraction of sig odor resps that are picked up from perfect cooperativity and re-distributed randomly
                                                %for -1 to 0; defines the fraction of sig odor resps that are picked up from a systematically distributed set of responses and re-distributed randomly
%decoder parameters
duration = 1;             %odor stim duration number to use for analysis (1 - 1s, 2 - 20s, 3 - 60s)
integration_window = 3;  %in s, the duration from stim_onset over which the dF/F traces are averaged to calculate response size 

n_sim_reps = 30;
%running simulation repeatedly to 
var_vec = [.1:.1:.8, .85, .9, .95];
score_mat = zeros(length(var_vec), n_sim_reps) + nan;
for var_n = 1:length(var_vec)
    m_sparseness = var_vec(var_n);
    
    for sim_rep_n = 1:n_sim_reps
        
        %function to generate re-sampled population responses
        [sim_data_mat, sp_vec] = resample_cell_pop(orig_data_mat{direc_n, 1}, n_cells, n_odors, m_sparseness, sd_sparseness, cooperativity);

        %function to take a re-sampled population response matrix and do a PCA and calculate odor separablity metrics
        classification_scores = odor_classifier_v2(sim_data_mat, duration, integration_window, sp_vec);
        c_score = classification_scores{1, 2};
        score_mat(var_n, sim_rep_n) = c_score;
    
        clear sim_data_mat
    end
    
end

