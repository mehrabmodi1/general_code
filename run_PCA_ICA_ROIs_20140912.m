%Written by Mehrab N Modi, 20140912. The purpose of this code is to call
%the functions in the PCA-ICA image segmentation toolbox associated with
%the paper Eran Mukamel, Axel Nimmerjahn and Mark Schnitzer, NEURON (2009).

clear all
close all


%House-keeping
%-----------------------

path('C:\Stuff\CSHL\Glenn lab\Code\Mukamel_PCA_ICA_ROIs\', path);       %adding folder with Mukamel code to Matlab's list of function paths
direc = 'C:\Data\CSHL\20140908\odor_stim_set2\';                        %directory with raw data in it
save_direc = 'C:\Data\CSHL\Analysed_data\';                             %main directory to save processed data files

[no_frames, no_trials, frame_rate, zoom, d_path, f_name, tr_tg_no] = tiff_info(direc);
frame_time = 1./frame_rate;                                             %in s

save_path = [save_direc d_path];                                        %final sub-directory to save processed data files
mkdir(save_path);

%generating mean fluorescence image
trial_no = 4;
trial_no_f = trial_no + tr_tg_no - 1;           %adding the trial tag number of the first trial of this set to get the right filename
trial_no_f = sprintf('%03.0f', trial_no_f);
filepath = [direc f_name trial_no_f, '.tif'];  %full file path for currently analysed trial
for frame_no = 1:30
    frame = imread(filepath, frame_no);
    frame = double(frame);
    if frame_no == 1
        ave_frame = zeros(size(frame));
    elseif frame_no > 1
        ave_frame = ave_frame + frame;
    else
    end
end
ave_frame = ave_frame./30;



%running Mukamel Code
%-----------------------

nPCs = 200;
dsamp = [1, 1];                 %downsampling factors for spatial and temporal sample rates
tem_weight = 0.5;               %no. between 0 and 1 specifying weight of temporal information in spatio-temporal ICA
termtol = 0.05;                 %fractional change in output at which to terminate iteration of 'fixed point algorithm'
maxrounds = 2000;               %max no of iterations to run
disp_mode = 'contour';          %'contour' shows contours of all spatial filters. 'series', each one at a time

%loop to read in multiple trials
%for trial_no = 1:3
    trial_no = 4;
    trial_no_f = trial_no + tr_tg_no - 1;           %adding the trial tag number of the first trial of this set to get the right filename
    trial_no_f = sprintf('%03.0f', trial_no_f);
    file_path = [direc f_name trial_no_f, '.tif'];  %full file path for currently analysed trial

    %step 1: CellsortPCA
    [mixedsig, mixedfilters, CovEvals, covtrace, movm,movtm] = CellsortPCA(file_path, [], nPCs, dsamp, save_path, []);
    
    %step 2: CellsortChoosePCs                  
    [PCuse] = CellsortChoosePCs(file_path, mixedfilters)
    nIC = floor(length(PCuse)./2);                       %no. of ICs to extract
    
    %step 3: CellsortPlotPCspectrum
    CellsortPlotPCspectrum(file_path, CovEvals, PCuse)
    
    %step 4: CellsortICA
    ica_A_guess = rand(length(PCuse), nIC);           %*********************
    [ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals, PCuse, tem_weight, nIC, ica_A_guess, termtol, maxrounds);
    
    %step 5: CellsortICAplot
    CellsortICAplot(disp_mode, ica_filters, ica_sig, ave_frame, [0, (frame_time.*no_frames)], frame_time, 1, 1, nIC)
    
%end

