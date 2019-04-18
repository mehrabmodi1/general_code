clear all
close all

list_direc = ['C:\Data\CSHL\dataset_list_Toshi_KC_DN_Led_20150708.txt'];
fid = fopen(list_direc);
direc_counter = 0;

color_vec = [[.35, .55, .75]; [.75, .55, .35]; [.55, .35, .75]; [.35, .55, .75]; [.75, .55, .35]; [.75, .35, .55]];
tr_list = [3:10, 13:20];


%figure property initialisation variables
plot_height = 300;
plot_width = 300;
axis_font_size = 15;
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
        
    %loading extracted raw fluorescence data matrices written by
    %raw_dff_extractor
    dataset = load([direc 'expt.mat']);
    dataset = dataset.data;
    stim_time = dataset(1).stim.stimLatency.*1000;              %stimulus onset time in ms
    stim_time = stim_time + 625;                                %added delay from valve opening to odor at pipe outlet
    frame_time = dataset(1).info.framePeriod .* 1000;           %frame time in ms
    stim_frame = floor(stim_time./frame_time);                  %frame no at which odor reached fly
    
    raw_data_mat = load([direc 'expt_raw_traces.mat']);
    raw_data_mat = raw_data_mat.raw_data_mat;
    
   
    %calculating dF/F traces and creating the sparse, 4-D, nan-filled
    %dff_data_mat 
    [dff_data_mat, stim_mat, prot_switch_trials] = cal_dff_traces_20150710(raw_data_mat, dataset, list_direc);
    clear raw_data_mat
    
    %building list of dropped trials
    nan_list = isnan(squeeze(dff_data_mat(1, 1, :, :)));
    nan_list_compressed = sum(nan_list, 2);
    n_odors = size(dff_data_mat, 4);
    dropped_trials = find(nan_list_compressed == n_odors);
    
    tr_ni = 0;
    while tr_ni < length(tr_list)
        tr_ni = tr_ni + 1;
        trial_n = tr_list(tr_ni);
        
        %skipping trial if its been dropped from analysis
        if isempty(intersect(trial_n, dropped_trials)) == 0
            continue
        else
        end
        
        curr_stack = dataset(trial_n).imageStack;
        
        
        baseline_frame = nanmean(curr_stack(:, :, (stim_frame - 22):(stim_frame - 2)), 3);
        baseline_stack = repmat(baseline_frame, [1, 1, size(curr_stack, 3)]);
        dff_stack = (curr_stack - baseline_stack)./baseline_stack;
        
        mean_dff_frame = mean(dff_stack, 3);
        figure(1)
        try
            overlayDFFandBaseline(dataset(trial_n))
        catch
            keyboard
        end
        disp(['trial number ', int2str(trial_n)])
        a = input('2 to skip to post-pairing trials, 1 for more trials, 0 for next dataset')
        
        if a == 0
            close figure 1
            break
        else
        end
        if a == 2
            tr_ni = 9;
        else
        end
        
        close figure 1
    end

end
fclose(fid);