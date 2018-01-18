clear all
close all

direc_lists_mat = [...
                {'C:\Data\Data\Raw_data\dataset_lists\dataset_list_axon_train_stim_20180118.xls'}...
                    ];
                
n_direc_lists = size(direc_lists_mat, 1);

%loop to go through all directory lists
for direc_list_n = 1:n_direc_lists
    list_direc = direc_lists_mat{direc_list_n, 1};
    [del, curr_direc_list] = xlsread(list_direc, 1);
    n_dirs = size(curr_direc_list, 1);
    direc_counter = 0;
    
    %parsing direc list path for name of direc list
    namei = findstr(list_direc, '\');
    namei = namei(end) + 1;
    dir_list_name = (list_direc(namei:(end - 4)));
    
    %loop to go through all experiment datasets listed in list file
    for direc_counter = 1:n_dirs
        %% House-keeping
        direc = curr_direc_list{direc_counter, 1};
        direc = [direc, '\'];
        prev_direc = pwd;
        cd([direc]);
        tif_list = dir('*.tif');
        
        ref_im = ScanImageTiffReader([direc, tif_list(1).name]).data();
        ref_im = mean(ref_im, 3)';
        
        %loading in manually drawn, FIJI ROIs
        if exist([direc, 'ROIs']) ~= 7
            mkdir([direc, 'ROIs']);
        	unzip([direc, 'RoiSet.zip'], [direc, 'ROIs\']);     %extracting ROIs from .zip file
        else
        end
        prev_direc = pwd;
        cd([direc, 'ROIs\']);
        ROI_list = dir('*.roi');
        cd(prev_direc);
        n_ROIs = size(ROI_list, 1);
        
        ROI_mat = zeros(size(ref_im, 1), size(ref_im, 2), n_ROIs);
        for ROI_n = 1:n_ROIs
            curr_name = ROI_list(ROI_n).name;
            curr_ROI = ReadImageJROI([direc, 'ROIs\', curr_name]);
            [test_pointsx, test_pointsy] = find(ROI_mat(:, :, ROI_n) == 0);
            [in_pts, on_pts]= inpolygon(test_pointsx, test_pointsy, curr_ROI.mnCoordinates(:, 1), curr_ROI.mnCoordinates(:, 2));
            in_points = find(in_pts == 1);
            in_points = [in_points; find(on_pts == 1)];
            ROI_mat(test_pointsy(in_points), test_pointsx(in_points), ROI_n) = 1;
        end
        keyboard
        
    end
end