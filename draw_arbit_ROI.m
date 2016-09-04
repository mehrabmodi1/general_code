clear all
close all

saving = 1;
fid=fopen('C:\Data\CSHL\dataset_list_20141208.txt');
while 1
    
    direc = fgetl(fid);
    if ~ischar(direc),   break,   end
    
    direc = [direc '\'];
    %House-keeping
    %-----------------------
    save_direc = 'C:\Data\CSHL\Analysed_data\';                             %main directory to save processed data files

    [no_frames, no_trials, frame_rate, zoom, d_path, f_name, tr_tg_no, ave_frame] = tiff_info(direc);
    [isi, odor_t_list, odorNames, n_reps, stim_time] = load_stim_param(direc);

    frame_time = 1./frame_rate;                                             %in s
    stim_frame = floor(stim_time(1, 1)./frame_time);                              %frame number at which stim was delivered

    save_path = [save_direc d_path];                                        %final sub-directory to save processed data files
    mkdir(save_path);

    trial_no = 1;
    trial_no_f = trial_no + tr_tg_no - 1;           %adding the trial tag number of the first trial of this set to get the right filename
    trial_no_f = sprintf('%03.0f', trial_no_f);
    file_path = [direc f_name trial_no_f, '.tif'];  %full file path for currently analysed trial


    %reading in trial 1 frames
    for frame_n = 1:30
        frame = double(imread(file_path, frame_n));

        if frame_n == 1
            ave_frame = frame;
        else
        end

        if frame_n > 1
            ave_frame = ave_frame + frame;
        else
        end

    end

    ave_frame = ave_frame./30;

    colormap('gray')
    imagesc(ave_frame)
    ROI = roipoly;

    ROI_int_mat = zeros(size(ROI));
    ROI_int_mat = ROI_int_mat + ROI;
    
    %saving
    if saving == 1
        save([save_path 'ROI_int_mat.txt'], 'ROI_int_mat', '-ASCII');
        save([direc 'ROI_mat.txt'], 'ROI_int_mat', '-ASCII');
    else
    end
    
    disp(direc)
end

fclose(fid);