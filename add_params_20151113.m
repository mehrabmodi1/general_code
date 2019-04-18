clear all
close all

list_direc = ['D:\Data\CSHL\dataset_list_dopamine_20150421.txt']; %expt datasets

blocks2compare = [1, 2];            %Specify the two blocks between which to compare KC pop responses
treatment_odor = 2;                 %The odor paired with something that is to be compared with other odors

an_trial_window = 5;                %Specify the window of trials to use from each block, for each odor (eg. first 4 trials)    
fid = fopen(list_direc);
direc_counter = 0;

plot_summary_an_only = 1;       %1 - only plots analyses across datasets, 0 - also plots for individual datasets
plot_PC_trajs = 1;              %if plot_summary_an_only is 0, this comes into effect
summary_statistic = 1;          % 0 - use mean; 1 - use median

sig_cell_block_mat = [];

all_c_dists = [];
mean_c_dists_all = [];
c_dist_lengths = [];
all_block1_dists = [];
color_vec = [[.75, .35, .25]; [.50, .70, .30]; [.55, .35, .75]; [.35, .55, .75]; [.75, .55, .35]; [.75, .35, .55]];
n_cells_saved = [];
mean_vec_first_pooled = [];
mean_vec_last_pooled = [];
se_vec_first_pooled = [];
se_vec_last_pooled = [];
odor1_corrs = [];
odor2_corrs = [];
all_c_dists_cosine = [];

%figure property initialisation variables
plot_height = 200;
plot_width = 200;
axis_font_size = 8;
tick_length = [0.02, 0.05];
line_widths = 2;
marker_sizes = 1;
marker_sizes_f = 4;
box_l_width = 0.5;

centroid_dists_saved = [];


%loop to go through all experiment datasets listed in list file
while 1
    direc_counter = direc_counter + 1;
    direc = fgetl(fid);
       
    if ischar(direc) ~= 1
        break
    else
    end
    
    direc = [direc '\'];
        
    %replacing C: with D:
    a = findstr(direc, 'C:');
    direc(a) = 'D';

    
    
    %loading extracted raw fluorescence data matrices written by
    %raw_dff_extractor
    dataset = load([direc 'expt.mat']);
    dataset = dataset.data;
    keyboard
    
    
    %saving new, combined mat file and continuing analysis with it
    data = dataset;
    save([direc 'expt'], 'data');
    clear data
    cd(direc)

    
end
fclose(fid);






