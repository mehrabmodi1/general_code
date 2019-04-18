clear all
close all

cell_direc_list_path = 'C:\Data\Data\Raw_data\dataset_lists\dataset_list_sim_patch_2P_cells_nlsGC6f.xls';

[del, cell_direc_list] = xlsread(cell_direc_list_path, 1);
n_cells = size(cell_direc_list, 1);

for cell_n = 1:n_cells
    curr_cell_direc = cell_direc_list{cell_n, 1}; 
    prev_direc = pwd;
    cd(curr_cell_direc);
    dir_contents_h5 = dir_date_sorted(curr_cell_direc, '*.h5');
    dir_contents_tif = dir_date_sorted(curr_cell_direc, '*.tif');
    Vscale = 10;
    
    %analysing current-injection trials
    %% Reading in and aligning current clamp traces
    patch_fname_buzz = arrayfun(@(a) not(isempty(strfind(a.name, 'buzz'))), dir_contents_h5);   %Generates a list of indices of dir_contents_h5.name that contain the string 'buzz'
    num_vec = 1:1:size(patch_fname_buzz, 2);
    patch_fname_buzz = num_vec(patch_fname_buzz);               %converted vector of logicals to vector of trial numbers
    if length(patch_fname_buzz) > 1
        fname_final = 1;
        patch_fname_final =  arrayfun(@(a) not(isempty(strfind(a.name, 'final'))), dir_contents_h5);
        num_vec = 1:1:size(patch_fname_final, 2);
        patch_fname_final = num_vec(patch_fname_final);
    
        patch_fname_buzz = intersect(patch_fname_buzz, patch_fname_final);
    else
        fname_final = 0;
    end
    
    try
        raw_traces = ws.loadDataFile(dir_contents_h5(patch_fname_buzz).name);
    catch
        keyboard
    end
    
    n_trials = raw_traces.header.NSweepsPerRun;
    tr_duation = raw_traces.header.SweepDuration;
    sample_rate = raw_traces.header.AcquisitionSampleRate;
    
    raw_trace_mat = zeros((tr_duation.*sample_rate), n_trials) + nan;
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
    n_trials = n_trials - length(skipped_tr_list);
    
    %aligning traces in time
    thresh = -3.0;
    al_trace_mat = zeros(6.*sample_rate, n_trials);
    saved_onsets = zeros(n_trials, 1);
    for tr_n = 1:n_trials
        curr_tr = raw_trace_mat(:, tr_n);
        buzz_onseti = find(curr_tr > thresh);
        buzz_onseti = buzz_onseti(1);
        saved_onsets(tr_n) = buzz_onseti./sample_rate;
       
        al_trace_mat(:, tr_n) = curr_tr(buzz_onseti - (3.*sample_rate):(buzz_onseti + (3.*sample_rate) - 1)).*Vscale;
        
    end
    
    %% Reading in fluorescence traces
    %loading first .tiff and asking for a manual ROI
    twop_fname_buzz = arrayfun(@(a) not(isempty(strfind(a.name, 'buzz'))), dir_contents_tif);   %Generates a list of indices of dir_contents_tif.name that contain the string 'buzz'
    num_vec = 1:1:size(twop_fname_buzz, 2);
    twop_fname_buzz = num_vec(twop_fname_buzz);             %converted vector of logicals to vector of trial numbers
    
    if fname_final == 1
        twop_fname_final =  arrayfun(@(a) not(isempty(strfind(a.name, 'final'))), dir_contents_tif);
        num_vec = 1:1:size(twop_fname_final, 2);
        twop_fname_final = num_vec(twop_fname_final);
            
        twop_fname_buzz = intersect(twop_fname_buzz, twop_fname_final);
    else
    end
    
    raw_F_mat = [];
    for trial_n = 1:n_trials
        stack_obj = ScanImageTiffReader([curr_cell_direc, '\', dir_contents_tif(twop_fname_buzz(trial_n)).name]);
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
    
    %aligning F traces to i-clamp traces
    al_F_mat = zeros(ceil(6./frame_time), n_trials);
    for trial_n = 1:n_trials
        onset_fr = ceil(saved_onsets(trial_n)./frame_time);
        n_frs_2s = ceil(3./frame_time); 
        al_F_mat(:, trial_n) = dff_mat( (onset_fr - n_frs_2s + 1):(onset_fr + n_frs_2s - 1), trial_n);
    end
    
    
    %% Plotting
    %aligned i-clamp traces
    figure(1)
    margin_val_h = 0.12;
    margin_val_v = 0.10;
    subplot_tight(2, 1, 1, [margin_val_h, margin_val_v])
    t_vec  = (0:(1./sample_rate):6)';
    t_vec = t_vec - 3;
    t_vec(end) = [];
    ave_trace = mean(al_trace_mat, 2);
    plot(repmat(t_vec, 1, n_trials), al_trace_mat, 'Color', [0.6, 0.6, 0.6])
    hold on
    plot(t_vec, ave_trace, 'Color', [0, 0, 0], 'LineWidth', 2)
    xlabel('time (s)')
    ylabel('membrane voltage (mV)')
    hold off
    
    subplot_tight(2, 1, 2, [margin_val_h, margin_val_v])
    t_vec  = (0:frame_time:6)';
    t_vec = t_vec - 3;
    ave_trace = mean(al_F_mat, 2);
    plot(repmat(t_vec, 1, n_trials), al_F_mat, 'Color', [0.6, 0.6, 0.6])
    hold on
    plot(t_vec, ave_trace, 'Color', [0, 0, 0], 'LineWidth', 2)
    xlabel('time (s)')
    ylabel('dF/F')
    hold off
    
   keyboard
   try
    close figure 1
   catch
   end
        
end