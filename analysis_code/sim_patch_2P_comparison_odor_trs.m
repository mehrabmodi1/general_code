clear all
close all

color_vec = load('C:\Data\Code\general_code\std_color_vec.txt');
cell_direc_list_path = 'C:\Data\Data\Raw_data\dataset_lists\dataset_list_sim_patch_2P_odor_resps.xls';

[del, cell_direc_list] = xlsread(cell_direc_list_path, 1);
n_cells = size(cell_direc_list, 1);

for cell_n = 1:n_cells
    curr_cell_direc = cell_direc_list{cell_n, 1}; 
    curr_cell_direc = [curr_cell_direc, '\'];
    prev_direc = pwd;
    cd(curr_cell_direc);
    remove_small_tifs(curr_cell_direc);
    [stim_mat, stim_mat_simple, column_heads, color_vec] = load_params_trains(curr_cell_direc, []);
    
    dir_contents_h5 = dir_date_sorted(curr_cell_direc, '*.h5');
    dir_contents_tif = dir_date_sorted(curr_cell_direc, '*.tif');
    Vscale = 10;
    
    if exist([curr_cell_direc, 'al_trace_mat.mat' ]) ~= 2
        %analysing current-injection trials
        %% Reading in and aligning current clamp traces
        try
            raw_traces = ws.loadDataFile(dir_contents_h5(1).name);
        catch
            keyboard
        end

        n_trials = length(fieldnames(raw_traces)) - 1;
        tr_duration = raw_traces.header.SweepDuration;
        sample_rate = raw_traces.header.AcquisitionSampleRate;

        raw_trace_mat = zeros((tr_duration.*sample_rate), n_trials) + nan;
        skipped_tr_list = [];
        for tr_n = 1:n_trials
            str_trn = num2str(tr_n);
            pad = zeros(1, (4 - size(str_trn, 2)));
            padded_n = [num2str(pad), str_trn];
            spacesi = findstr(padded_n, ' ');
            padded_n(spacesi) = [];
            try
                eval(['raw_trace_mat(:, tr_n) = raw_traces.sweep_', padded_n, '.analogScans(:, 2);']);
            catch
                skipped_tr_list  = [skipped_tr_list; tr_n];
            end
        end
        raw_trace_mat(:, skipped_tr_list) = [];
        n_trials = size(raw_trace_mat, 2);

        %sampling traces around odor delivery period
        al_trace_mat = zeros(25.*sample_rate, n_trials);
        for tr_n = 1:n_trials
            stim_sample = floor(stim_mat_simple(tr_n, 7).*sample_rate);
            curr_tr = raw_trace_mat(:, tr_n);
            al_trace_mat(:, tr_n) = curr_tr(stim_sample - (5.*sample_rate):(stim_sample + (20.*sample_rate) - 1)).*Vscale;

        end

        %% Reading in fluorescence traces
        %loading first .tiff and asking for a manual ROI
        raw_F_mat = [];
        for trial_n = 1:n_trials
            stack_obj = ScanImageTiffReader([curr_cell_direc, '\', dir_contents_tif(trial_n).name]);
            stack = stack_obj.data();
            stack = permute(stack,[2 1 3]);
            n_frames = size(stack, 3);
            [frame_time, zoom, n_chans] = SI_tif_info(stack_obj);
            if n_chans == 2
                g_frames = 1:2:(n_frames - 1);
                stack = stack(:, :, g_frames);
                n_frames = size(stack, 3);
           
            else
            end



            %getting ROI
            if trial_n == 1
                if exist([curr_cell_direc, '\', 'ROI.mat'], 'file') == 2
                    ROI_mat = load([curr_cell_direc, '\', 'ROI.mat']);
                    ROI_mat = ROI_mat.ROI_mat;
                else
                    ave_im = mean(stack, 3);
                    figure(1)
                    imagesc(ave_im)
                    ROI_mat = roipoly();
                    close figure 1
                    save([curr_cell_direc, '\', 'ROI.mat'], 'ROI_mat');
                end
            else
            end
            curr_pixi = find(ROI_mat == 1);


            %extracting fluorescence traces
            F_trace = zeros(n_frames, 1);
            for frame_n = 1:n_frames
                curr_fr = stack(:, :, frame_n);
                F_trace(frame_n, 1) = mean(curr_fr(curr_pixi));
            end
            raw_F_mat = pad_n_concatenate(raw_F_mat, F_trace, 2, nan);

        end
        %calculating dF/F
        baselines = mean(raw_F_mat(1:(ceil(3./frame_time)), :), 1);
        baselines = repmat(baselines, size(raw_F_mat, 1), 1);
        dff_mat = (raw_F_mat - baselines)./baselines;

        %sampling F traces around odor period
        al_F_mat = zeros(ceil(25./frame_time), n_trials);

        for trial_n = 1:n_trials
            stim_fr = floor(stim_mat_simple(trial_n, 7)./frame_time);
            n_frs_5s = ceil(5./frame_time); 
            n_frs_20s = ceil(20./frame_time);
            al_F_mat(:, trial_n) = dff_mat( (stim_fr - n_frs_5s + 1):(stim_fr + n_frs_20s - 1), trial_n);
        end

        al_trace_mat_full = al_trace_mat;
        al_F_mat_full = al_F_mat;

        %saving data to file
        save([curr_cell_direc, 'al_trace_mat.mat' ], 'al_trace_mat_full');
        save([curr_cell_direc, 'al_F_mat.mat' ], 'al_F_mat_full');
    elseif exist([curr_cell_direc, 'al_trace_mat.mat' ]) == 2
        al_trace_mat_full = load([curr_cell_direc, 'al_trace_mat.mat' ]);
        al_trace_mat_full = al_trace_mat_full.al_trace_mat_full;
        al_F_mat_full = load([curr_cell_direc, 'al_F_mat.mat' ]);  
        al_F_mat_full = al_F_mat_full.al_F_mat_full;
        
        %reading in frame time from the first trial
        stack_obj = ScanImageTiffReader([curr_cell_direc, '\', dir_contents_tif(1).name]);
        [frame_time, zoom, n_chans] = SI_tif_info(stack_obj);
        
        raw_traces = ws.loadDataFile(dir_contents_h5(1).name);
        sample_rate = raw_traces.header.AcquisitionSampleRate;
        stim_sample = floor(stim_mat_simple(1, 7).*sample_rate);
       
        n_trials_all = size(al_F_mat_full, 2);
    end
    
    %% Plotting
    odor_list = unique(stim_mat_simple(:, 2));
    for odor_n = 1:length(odor_list)
        odor_ni = odor_list(odor_n);
        curr_trs = find(stim_mat_simple(:, 2) == odor_ni);
        curr_trs(curr_trs > n_trials_all) = [];
        al_trace_mat = al_trace_mat_full(:, curr_trs);
        al_F_mat = al_F_mat_full(:, curr_trs);
        n_trials = length(curr_trs);
                
        stim_fr = floor(stim_mat_simple(curr_trs(1), 7)./frame_time);
        n_frs_5s = ceil(5./frame_time); 
        n_frs_20s = ceil(20./frame_time);
        stim_frs = compute_stim_frs(stim_mat, curr_trs(1), frame_time);
        stim_frs = stim_frs - (stim_fr - n_frs_5s + 1);
        stim_samples = compute_stim_frs(stim_mat, curr_trs(1), 1./sample_rate);
        stim_samples = stim_samples - (stim_sample - (5.*sample_rate) );
        
        
        %aligned i-clamp traces
        figure(1)
        ave_trace = mean(al_trace_mat, 2);
        plot(al_trace_mat, 'Color', [0.6, 0.6, 0.6])
        hold on
        plot(ave_trace, 'Color', [0, 0, 0], 'LineWidth', 2)
        xlabel('time (s)')
        ylabel('membrane voltage (mV)')
        set_xlabels_time(1, (1./sample_rate), 5)
        if odor_ni <= size(color_vec, 1)
            color_n = odor_ni;
        else
            color_n = rem(odor_ni, size(color_vec, 1) );
        end
        add_stim_bar(1, stim_samples, color_vec(color_n, :))
        
        hold off
        
        figure(2)
        ave_trace = mean(al_F_mat, 2);
        plot(al_F_mat, 'Color', [0.6, 0.6, 0.6])
        hold on
        plot(ave_trace, 'Color', [0, 0, 0], 'LineWidth', 2)
        xlabel('time (s)')
        ylabel('dF/F')
        set_xlabels_time(2, frame_time,5)
        add_stim_bar(2, stim_frs, color_vec(color_n, :))
        
        hold off

       %keyboard
       del = input('');
       try
        close figure 1
        close figure 2
       catch
       end
       
    end
   
        
end