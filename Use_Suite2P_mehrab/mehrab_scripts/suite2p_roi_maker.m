clear all
close all

results_direc = 'E:\Data\Analysed_data\Suite2p\Results\';
results_direc_manual_ROIs = 'E:\Data\Analysed_data\Manual_ROIs\';
reg_tif_direc = 'E:\Data\Analysed_data\Suite2p\Reg_Tiff\';
raw_direc_base = 'E:\Data\Raw_Data_Current\Resonant\';

extract_both_channels = 0;

%% Reading in manually created direc list
direc_list = 'E:\Data\Raw_Data_Current\dataset_lists\KC_d5HT1b_PABAEL_201908set.xls';
[del, raw_direc_list] = xlsread(direc_list, 1);
raw_direc_list_copy = raw_direc_list;

%% Setting up loop to repeatedly run Suite2P for each identified dataset directory (ie. all sub-directories of raw_direc_base that contain stim params files.)
n_direcs_analysed = 0;
rem_dir_list = [];
for raw_direc_n = 1:size(raw_direc_list, 1)
   
    % House-keeping
    direc = raw_direc_list{raw_direc_n, 1};
    direc = [direc, '\'];
    remove_small_tifs(direc);
    prev_direc = pwd;
    cd([direc]);
    dataset_namei = findstr(direc, '\20');
    raw_direc = direc((dataset_namei + 1):end);
    raw_direc_list{raw_direc_n, 1} = raw_direc;    
    %% Pre-prep
    %raw_direc = [raw_direc '\'];
    
    %Checking if this directory has already been analysed with Suite2P or ManualROIextracter

    if isdir([results_direc_manual_ROIs, raw_direc]) == 1
        %raw_direc_list(raw_direc_n) = [];
        continue
    else
    end
    
    if isdir([results_direc raw_direc]) == 0
        prev_direc = pwd;
        cd([raw_direc_base, raw_direc]);
        disp(['currently analysing ' raw_direc_base, raw_direc]);
        
        %Getting rid of empty .tiff files in raw data folder(s) left behind by ScanImage
        remove_small_tifs([raw_direc_base, raw_direc]);

        %setting up directory structure for Suite2P
        raw_direc = setup_raw_direc(raw_direc_base, raw_direc);
        raw_direc_list{raw_direc_n, 1} = raw_direc;
        
        %skipping if current directory doesn't meet criteria in setup_raw_direc
        if isempty(raw_direc) == 1
            continue
        else
        end
        
        %setting up master_file.m and make_db.m
        ft_factor = setup_Suite2P_files(raw_direc_base, raw_direc, 'E:\Code\general_code_repo2\general-code\Use_Suite2P_mehrab\mehrab_scripts\Suite2P_starter_files\master_file.m');
        
        prev_dir = pwd;
        cd([raw_direc_base, raw_direc]);
        %running Suite2P
%         try
            master_file
%         catch
%             keyboard
%         end
        n_direcs_analysed = n_direcs_analysed + 1;

        cd(prev_direc);
    else
        disp([raw_direc_base, raw_direc ' has already been analysed... skipping.'])
        rem_dir_list = [rem_dir_list; raw_direc_n];
    end
    disp(['analysed ' int2str(n_direcs_analysed) ' new datasets.']);
end
%raw_direc_list(rem_dir_list) = [];



%% Loop to repeatedly call ROI_prune to manually curate ROIs identified by Suite2P
for raw_direc_n = 1:size(raw_direc_list, 1)
    raw_direc = raw_direc_list{raw_direc_n, 1};
    raw_direc = raw_direc_with_1(raw_direc_base, raw_direc);
    disp(raw_direc);

    %% Loading in Suite2P analysis results file
    prev_direc = pwd;
    cd([results_direc, raw_direc])
    
    if exist([results_direc, raw_direc, '\ROIs_pruned.txt']) ~= 2
        new_main       %this is the ROI pruning GUI that comes with Suite2P
        keyboard       %don't remove this keyboard - it's needed for normal running.
        del = [];
        save([results_direc, raw_direc, '\ROIs_pruned.txt'], 'del');
    else
        disp([results_direc, raw_direc, ' has ROIs already pruned. Skipping.']);
    end
    
    %saving pruned Suite2P ROIs as a simple ROI.mat file.
    newest_results_file = find_newest_file([results_direc, raw_direc], '_proc');
    if isempty(newest_results_file) == 1
        continue
    else
    end
    disp(['Loading Suite2P results file ' newest_results_file])
    data_mat = load([results_direc, raw_direc, '\', newest_results_file]);

    try
        data_mat = data_mat.dat;
        disp('Done loading.')
    catch
        del = [];
        save([raw_direc_base, raw_direc, 'skip_direc.txt'], 'del');
        disp(['no data in ' results_direc, raw_direc, dir_contents(max_datenum(2)).name, '... skipping.']);
        continue
    end
    cd(prev_direc);

    %Creating ROI matrix that includes only ROIs manually classified as real cells
    [ROI_mat] = setup_Suite2P_ROIs(data_mat);
    save([results_direc, raw_direc, '\ROI_mat.mat'], 'ROI_mat');
    clear data_mat
    
    %checking if directory structure was extended with a \1\ folder at the
    %end and copying over the results file from ROI_prune to that folder
    copy_files_to_1_direc(results_direc, raw_direc);
    
end