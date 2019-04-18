%This program compares neuronal population responses across trials for the
%same odor stimulus, and across different odor stimuli.
%Mehrab Modi, 20141014

clear all
close all

list_path = 'C:\Data\CSHL\dataset_list_20141031.txt';

resp_frac_mat = [];

fid=fopen(list_path);
counter = 0;
while 1
    counter = counter + 1;
    direc = fgetl(fid);
    if ~ischar(direc),   break,   end
    direc = [direc '\'];
    
    
    disp(['dataset number ' int2str(counter)])


    %direc = 'C:\Data\CSHL\20141030\water_odor_stim\';

    %House-Keeping
    %-------------------------------------------------
    path(path, 'C:\Stuff\CSHL\Glenn lab\Code\')
    save_direc = 'C:\Data\CSHL\Analysed_data\';                             %main directory to save processed data files

    %reading in dataset parameters
    [n_frames, n_trials, frame_rate, zoom, d_path, f_name, tr_tg_num, ave_frame] = tiff_info(direc);
    [isi, odor_t_list, odorNames, n_reps, stim_time] = load_stim_param(direc);

    save_path = [save_direc d_path];                                        %final sub-directory to save processed data files
    mkdir(save_path);

    frame_time = 1./frame_rate;
    stim_frame = floor(stim_time(1, 1)./frame_time);
    odor_t_list = odor_t_list(1:n_trials);
    n_odors = max(odor_t_list);

    %reading in extracted raw data values
    raw_data_mati = load([save_path 'raw_data_mat.mat']);
    raw_data_mat = raw_data_mati.raw_data_mat;
    done_trials = raw_data_mat{1, 2};               %trial number until which data has been extracted already
    raw_data_mat = raw_data_mat{1, 1};

    n_cells = size(raw_data_mat, 2);

    %calculating baselines for each cell, followed by dF/F for each
    %fluorescence intensity value
    Fo_mat = mean(raw_data_mat(1:(stim_frame - 1), :, :), 1);
    Fo_mat = repmat(Fo_mat, [n_frames, 1, 1]);

    dff_data_mat = (raw_data_mat - Fo_mat)./Fo_mat;         %calculating dF/F

    %------------------------

    %calculating response sizes for each trial, for each cell and turning
    %each time-series into a single number
    response_mat = zeros(n_cells, n_trials) + nan;
    response_t_mat = zeros(n_cells, n_trials) + nan;
    p_val_mat = zeros(n_cells, n_trials) + nan;           %matrix of p values of ttest comparisons between response and baseline points, trial by trial
    for cell_n = 1:n_cells
        for trial_n = 1:n_trials
            response_vec = dff_data_mat(:, cell_n, trial_n);
            [pk_area, pk_fr, p] = find_pk_area(response_vec, stim_frame, frame_time); 
            response_mat(cell_n, trial_n) = pk_area;
            response_t_mat(cell_n, trial_n) = pk_fr.*frame_time;
            p_val_mat(cell_n, trial_n) = p; 
        end

    end

    sig_response_mat = abs(binarise_im(p_val_mat, 0.05) - 1);       %this is the matrix of significant response trials, 1's significant and 0's not.

    %--------------------------------

    %sorting response_mat, response_time_mat and sig_response_mat by odors
    sig_mat_sorted = sortrows([odor_t_list, ([1:1:n_trials])', sig_response_mat']);
    resp_mat_sorted = sortrows([odor_t_list, ([1:1:n_trials])', response_mat']);
    resp_t_sorted = sortrows([odor_t_list, ([1:1:n_trials])', response_t_mat']);

    sig_mat_sorted(:, 2) = [];
    resp_mat_sorted(:, 2) = [];
    resp_t_sorted(:, 2) = [];

    a = find(sig_mat_sorted == 0);
    resp_mat_sorted(a) = 0;
    resp_t_sorted(a) = 0;

    figure(1)
    imagesc(resp_mat_sorted')
    title('Areas under dF/F curve, normalised')

    figure(2)
    imagesc(sig_mat_sorted')
    title('t-test p-values')

    figure(3)
    imagesc(resp_t_sorted')
    title('Timing of peak in dF/F curve')

    %----------------------------------------------
    %classifying and counting responsive cells

    %reading in total number of cells counted manually from images
    n_cells_total = load([save_path 'n_cells_counted.txt']);


    %counting numbers of responsive cells to each odor
    odor_resp_vec = zeros(1, 4) + nan;
    for odor_n = 1:4
        a = find(odor_t_list == odor_n);
        if isempty(a) == 1
            continue
        else
        end

        curr_odor_trials = sig_response_mat(:, a);
        resp_t_vec = sum(curr_odor_trials, 2);
        resp_t_vec = resp_t_vec./length(a);
        responders = find(resp_t_vec > 0.5);

        frac_responders = length(responders)./n_cells_total;

        odor_resp_vec(odor_n) = frac_responders;
    end

    save([save_path 'responder_fractions.txt'], 'odor_resp_vec', '-ASCII');
    
    resp_frac_mat = [resp_frac_mat; odor_resp_vec];
    
end
fclose(fid);