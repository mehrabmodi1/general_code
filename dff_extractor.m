function [] = dff_extractor(direc, verbose)
%This program requires an ROI file to have been previously created and
%placed in the data directory. It then aligns each frame with a reference
%frame for each trial and then compares the reference frames across trials
%and plots the corrcoef of these.
%Mehrab N. Modi
path(path, 'C:\Stuff\CSHL\Glenn lab\Code\')

%direc = 'C:\Data\CSHL\20141027\water_odor_stim\';
save_direc = 'C:\Data\CSHL\Analysed_data\';                             %main directory to save processed data files

%reading in dataset parameters
[n_frames, n_trials, frame_rate, zoom, d_path, f_name, tr_tg_num, ave_frame] = tiff_info(direc);
[isi, odor_t_list, odorNames, n_reps, stim_time] = load_stim_param(direc);

%n_frames = 20;          %REMOVE THIS LINE!

save_path = [save_direc d_path];                                        %final sub-directory to save processed data files
mkdir(save_path);

frame_time = 1./frame_rate;
stim_frame = floor(stim_time(1, 1)./frame_time);
odor_t_list = odor_t_list(1:n_trials);
n_odors = max(odor_t_list);

ROI_mat = load([direc 'ROI_mat.txt']);
ROI_mat = im2bw(ROI_mat, 0.5);
ROI_labels = bwlabel(ROI_mat, 4);
n_cells = max(max(ROI_labels));

%checking if raw intensity data has already been extracted and saved to file
if exist([save_path 'raw_data_mat.mat']) == 2
    raw_data_mati = load([save_path 'raw_data_mat.mat']);
    raw_data_mat = raw_data_mati.raw_data_mat;
    done_trials = raw_data_mat{1, 2};               %trial number until which data has been extracted already
    raw_data_mat = raw_data_mat{1, 1};
    
    raw_data_mati = load([save_path 'raw_data_mat_unreg.mat']);
    raw_data_mat_unreg = raw_data_mati.raw_data_mat_unreg;
    raw_data_mat_unreg = raw_data_mat_unreg{1, 1};
    clear raw_data_mati
else
    
    raw_data_mat = zeros(n_frames, n_cells, n_trials) + nan;
    raw_data_mat_unreg = zeros(n_frames, n_cells, n_trials) + nan;          %raw data extracted without 2d xcorr
    done_trials = 0;
    
end

%picking up the thread and continuing with thre trials left to be extracted
if done_trials < n_trials
    %loop for trials
    for trial_n = (done_trials + 1):n_trials
        %building filename string
        trial_no_f = trial_n + tr_tg_num - 1;           %adding the trial tag number of the first trial of this set to get the right filename
        trial_no_f = sprintf('%03.0f', trial_no_f);
        file_path = [direc f_name trial_no_f, '.tif'];  %full file path for currently analysed trial

        %reading in 5 frames to calculate average image
        for frame_n = 1:5
            frame = double(imread(file_path, frame_n) );
            if frame_n == 1
                template_fr = frame;
            else
            end
            if frame_n > 1
                template_fr = template_fr + frame;
            else
            end
        end
        template_fr = template_fr./5;
        
        %==================
        %aligning this trial's template frame to that of trial 1
        if trial_n == 1
            template_fr1 = template_fr;
            save([save_path 'template_fr_t1.txt'], 'template_fr1', '-ASCII');
        else
        end
        
        if trial_n > 1
            template_fr1 = load([save_path 'template_fr_t1.txt']);
        else
        end
        
        
        c = normxcorr2(template_fr1, template_fr);
        template_fr = xcorr2_aligner(template_fr, c);
        
        %==================
        
        
        
        
        %aligning each frame to template and then reading out cell data
        for frame_n = 1:n_frames
            %building filename string
            trial_no_f = trial_n + tr_tg_num - 1;           %adding the trial tag number of the first trial of this set to get the right filename
            trial_no_f = sprintf('%03.0f', trial_no_f);
            file_path = [direc f_name trial_no_f, '.tif'];  %full file path for currently analysed trial
            frame = double(imread(file_path, frame_n) );
            %calculating 2D cross correlation
            c = normxcorr2(template_fr, frame);
            frame_orig = frame;
            frame = xcorr2_aligner(frame, c);
            
            %keeping all frames for this trial in memory to identify bad frames later
            if frame_n == 1
                saved_frame_mat = zeros(size(frame, 1), size(frame, 2), n_frames) + nan;
            else
            end
            saved_frame_mat(:, :, frame_n) = frame;

            %reading in intensities of all pixels of each cell's ROI for this frame
            for cell_n = 1:n_cells
                curr_cell_pixi = find(ROI_labels == cell_n);
                curr_cell_pix = frame(curr_cell_pixi);
                raw_data_mat(frame_n, cell_n, trial_n) = nanmean(curr_cell_pix); 
                curr_cell_pix = frame_orig(curr_cell_pixi);
                raw_data_mat_unreg(frame_n, cell_n, trial_n) = nanmean(curr_cell_pix); 
            end
            
        end
        disp(['trials done: ' int2str(trial_n)])

        %=================
%         figure(1)
%         subplot(2, 1, 1)
%         colormap('gray')
%         imagesc(template_fr)
%         title('template of trial1')
% 
%         subplot(2, 1, 2)
%         colormap('gray')
%         imagesc(template_fr)
%         title(['template of trial ' int2str(trial_n)])
%         
%         
        
        %=================
        
        %identifying bad frames and replacing data for those frames with nans
        %bad_frames = bad_frame_finder(saved_frame_mat, template_fr);       no good way of doing this without a second imaging channel that has no activity 
        
        
        
        %saving extracted data
        raw_data_mat = {raw_data_mat, trial_n};
        save([save_path 'raw_data_mat.mat'], 'raw_data_mat')
        raw_data_mat = raw_data_mat{1, 1};

        raw_data_mat_unreg = {raw_data_mat_unreg, trial_n};
        save([save_path 'raw_data_mat_unreg.mat'], 'raw_data_mat_unreg')
        raw_data_mat_unreg = raw_data_mat_unreg{1, 1};
    
    end
    
    
    
    
else
end

if verbose == 1
    %calculating dF/F 
    %----------------

    %calculating baselines for each cell
    Fo_mat = mean(raw_data_mat(1:stim_frame, :, :), 1);
    Fo_mat = repmat(Fo_mat, [n_frames, 1, 1]);

    dff_data_mat = (raw_data_mat - Fo_mat)./Fo_mat;         %calculating dF/F


    Fo_mat_unreg = mean(raw_data_mat_unreg(1:stim_frame, :, :), 1);
    Fo_mat_unreg = repmat(Fo_mat_unreg, [n_frames, 1, 1]);

    dff_data_mat_unreg = (raw_data_mat_unreg - Fo_mat_unreg)./Fo_mat_unreg;         %calculating dF/F

    %end

    odor_list = {odorNames.odour};
    %displaying pop activity for each odor, looping through trials for that odor
    for odor_n = 1:n_odors
        t_list = find(odor_t_list == odor_n);
        curr_odor = odor_list{1, odor_n};
        for trial_n = 1:length(t_list)
            curr_tr = t_list(trial_n);

            figure(1)
            subplot(2, 1, 1)
            imagesc(dff_data_mat(:, :, curr_tr)');
            title([curr_odor ' repeat number ' int2str(curr_tr + tr_tg_num - 1) ])
            drawnow
            subplot(2, 1, 2)
            imagesc(dff_data_mat_unreg(:, :, curr_tr)');
            title([curr_odor ' repeat number ' int2str(curr_tr + tr_tg_num - 1) ])
            drawnow
            input('')
        end

    end
elseif verbose == 0
end
    
end