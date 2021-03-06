clear all
close all

direc_lists_mat = [...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_alpha1_lowUS.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_alpha1_lowUS_set2.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_alpha1_lowUS_set3.xls'};... 
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_gamma1_lowUS_set1.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_alpha1_lowUS_backward_ctrl.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Yoshi_THnull_G1pedc.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_gamma1_Chrim_stim_lifetime.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_gamma1_lowUS_MB085C.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_gamma1_lowUS_MB085C_epoxy.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\dataset_list_Yoshi_PaBaEl_MBON_DAN_gamma1_lowUS_MB085C_nansets'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_gamma1_lowUS_MB085C_new_set.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_gamma1_lowUS_MB085C_epoxy_short_session.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_gamma1_lowUS_MB085C_epoxy_short_session_noUS.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_gamma1_lowUS_MB085C_epoxy_short_session_noUS_shortCS.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_gamma1_lowUS_MB085C_epoxy_short_session_low_LED.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_gamma1_lowUS_MB085C_epoxy_short_session_low_LED_1Hz.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_gamma1_lowUS_MB085C_epoxy_short_session_low_LED_0.5Hz.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_gamma1_epoxy_short_session_lowerLED_0.1Hz.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\fore_distr_MBONAlpha1_set3_highLED.xls'};...                    
                    %{'E:\Data\Raw_Data_Current\dataset_lists\fore_distr_MBONGamma2.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBON_DAN_alpha1_eligibility_trace_Yoshi.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MPPC_KC_set_final.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MPPC_mouse_set.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\KC_d5HT1b_PABAEL_201908set.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\KC_c739_PABAEL_201908set.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\fore_distr_MBONGamma2_set1.xls'};...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\fore_distr_MBONGamma2_set2.xls'};...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\PABAEL_MBONGamma2_set1.xls'};...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\PABAEL_MBONGamma2_set2_0.1Hz.xls'};...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\fore_distr_MBONGamma2_set3_nodistr_ctrl.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\KC_dense_plasticity.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\KC_sparse_plasticity_set1.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\KC_sparse_odor_resps_set1.xls'};...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_starved_Ara_PaBaEl.xls'};...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover.xls'};...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\fore_distr_MBONAlpha2sc.xls'};...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_simple_starved.xls'};...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\Alpha1_60strace_71C03LxA_MB043CGal4.xls'};...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\Alpha1_60strace_72D01LxAChr88tdT.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_starved_set2.xls'}...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\KC_c739_PABAEL_201908set.xls'}...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_starved_halfAra.xls'}...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_starved36_halfAra.xls'}...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_starved_halfAra_prehabituated.xls'}...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_starved_halfAra_prehabituated_strongUS.xls'}...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\Alpha1_60strace_71C03LxA_MB043CGal4_noChrisoncontrol.xls'}...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_starved_halfAra_prehabituated_strongUS_EL_handover.xls'}...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_starved_halfAra_prehabituated_strongUS_EL_second.xls'}...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\c739KC_PaBaEl_handover_prehabituated.xls'}...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\Alpha1_60strace_71C03LxA_MB043CGal4_Chrison_noLED_control.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Alpha1_60strace_71C03LxA_MB043CGal4_with_Chrimson_LED.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBONA1_PaBaEl_handover_starved_halfAra_prehabituated_strongUS.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\c739KC_PaBaEl_handover_rtrain_prehabituated.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_starved_halfAra_prehabituated_noLED_control.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Alpha1_60strace_71C03LxA_MB043CGal4_with_Chrimson_LED.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Alpha1_10strace_71C03LxA_MB043CGal4_with_Chrimson_LED.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Alpha1_30sCS_60strace_71C03LxA_MB043CGal4_with_Chrimson_LED.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\c305aKC_PaBaEl_handover_rtrain_prehabituated.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\test_list.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_simp_pairing_Berry.xls'};...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_pairing_Berry.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_pairing_Berry_ELsecond.xls'};...
<<<<<<< Updated upstream
                    {'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_pairing_Berry_5sipi.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_pairing_Berry_15s_ipi.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_pairing_Berry_25sipi_higherLED.xls'};...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_pairing_Berry_longpulse2.xls'};...
                     %{'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_pairing_Berry_generalizatin.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_pairing_Berry_noLEDctrl.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_WTG2Ap1_onrig_fedshock.xls'};...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_WTG2Ap1_onrig_fedshock_30V.xls'};...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_WTG2Ap1_onrig_fedshock_10V.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_WTG2Ap1_onrig_fedshock_45V.xls'};...
                      
                    
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_onrig_shockcalib_13F02LexA_DA44.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_onrig_shock_13F02_DA44.xls'};...
                    
                    
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_MB298B_MBONG4-G1G2_GcaMP6f_starved.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_GH146_GCaMP6f_fed.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_GH146_lobes_fed.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_SS01240_DPM_starved.xls'};...
                    
                    
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_WTalpha3_fedshock_BAplus.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_WTalpha3_fedshock_PAplus.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_WTalpha3_fedshock_unpairedctrl.xls'};...
                    
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_WTG2Ap1_fedshock_BAplus_10xsp.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_WTG2Ap1_fedshock_BAplus_1x.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_WTG2Ap1_fedshock_unpairedctrl_1x.xls'};...
                    
                    
                    
                    
                     %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_13F02_gcaMP6f_starved.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBONA1_Berry_simple_trace0.xls'};...
                    
=======
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_pairing_Berry_15s_ipi.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\MBONG2_PaBaEl_handover_pairing_Berry_longpulse2.xls'};...
%                   {'E:\Data\Raw_Data_Current\dataset_lists\MBONA1_Berry_simple_trace0.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\Berry_handover_MB298B_MBONG4-G1G2_GcaMP6f_starved.xls'};...
                    
>>>>>>> Stashed changes
%                   {'E:\Data\Raw_Data_Current\dataset_lists\Alpha1_60strace_71C03LxA_MB043CGal4_arena_trained.xls'};...
%                   {'E:\Data\Raw_Data_Current\dataset_lists\c739KC_PaBaEl_handover_rtrain_prehabituated.xls'};...
%                   {'E:\Data\Raw_Data_Current\dataset_lists\KCd5HT1b_Gamma_PaBaEl_handover_rtrain_prehabituated.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\KC_AB_c739_vertlobe_longtrace.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\vGlut_MBON_GCaMP.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\13F02Gal4_x_ACh3_choline_sensor.xls'};...
%                     {'E:\Data\Raw_Data_Current\dataset_lists\Alpha1_0strace_71C03LxA_MB043CGal4_Berry.xls'};...
<<<<<<< Updated upstream
                    %{'E:\Data\Raw_Data_Current\dataset_lists\d5HT1b_similar_od_handovers.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\c305a_similar_od_handovers.xls'};...
                    %{'E:\Data\Raw_Data_Current\dataset_lists\c739_similar_od_handovers.xls'};...
=======
                    {'E:\Data\Raw_Data_Current\dataset_lists\d5HT1b_similar_od_handovers.xls'};...
>>>>>>> Stashed changes
                    %{'E:\Data\Raw_Data_Current\dataset_lists\ufGCaMPs_fast_imaging1.xls'};...
                        

                                            ];
                
n_direc_lists = size(direc_lists_mat, 1);

save_path_base_manual = 'E:\Data\Analysed_data\Manual_ROIs\';  
save_path_base_Suite2P = 'E:\Data\Analysed_data\Suite2p\Results\';


fuse_ROIs = 0;          %0-no, 1-yes. This flattens dim3 of ROI_mat in case there is a multi-patch neuron that needs to be considered as a single object.
dilate_ROIs = 0; %15;        %This is the number of pixels by which the manually drawn ROIs are dilated before data extraction.
remove_ROI_ovlap = 1;    %enabling this switch gets rid of an overlapping pixels between ROIs from all ROIs.   
warp_ROIs = 0;           %typically used for KC datasets - user-specified landmarks used to warp ROI matrix to align ROIs with cells as tissue distorts over session.
extract_both_channels = 0;  
extract_traces_only = 0; %Enabling this switch makes the program only extract traces from a .tif file using given ROIs. It doesn't look for stim params/PID files of any sort. useful for data acquired on other rigs. 
do_noisefit_subtraction = 1;    %Enabling this switch makes the pipeline fit slow, row-wise noise and subtract that from each frame.
do_bk_subtraction = 1;      %Enabling this switch makes the pipeline subtract an estimate of the background offset value computed from the manually drawn background ROI
force_re_extraction = 0;    %Enabling this switch forces the extraction program to re-extract data without re-doing the manual steps.
force_re_extr_manual = 0;   %Enabling this switch forces the user to re-do manual alignment, ROI adjustment steps followed by re-extraction of data.

if force_re_extraction == 1
    qstring = 'Are you sure? Force re-extraction of all dataset lists?';
    str1 = 'Yes';
    str2 = 'No';
    button = questdlg(qstring,' ',str1,str2, 'No');
    if strcmp(button, 'No') == 1
        force_re_extraction = 0;
    else
    end
else
end

%looping through all directory lists and all datasets once to save mean frames and again to save manually determined slow x-y  offsets for each trial, as well as determine bad z-drift trials
for do_over = 1:2
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
            direc1 = get_path_1_direc(direc);       %in case raw direc has \1 folders created by the Suite2P program
            %checking if ROI_type is manualROIs (0) or Suite2P (1)
            if length(direc) == length(direc1)
                ROI_type = 0 ;
            elseif length(direc) < length(direc1)
                ROI_type = 1;
                direc = direc1;
            else
            end
                    
            remove_small_tifs(direc);
           
            prev_direc = pwd;
            cd([direc]);
            dataset_namei = findstr(direc, '20');
            dataset_name = direc((dataset_namei):end);
            if ROI_type == 0
                save_path_base = save_path_base_manual;
            elseif ROI_type == 1
                save_path_base = save_path_base_Suite2P;
            else
            end
                
            save_path = [save_path_base, dataset_name, '\' ];
            save_path = get_path_1_direc(save_path);
           
            %checking if dataset has already been analysed
            old_dir = pwd;
            if isdir(save_path) == 1
                cd(save_path);

                tif_list = dir('*.tif');        %checking if trial 1 tif has been copied over to results folder
                if isempty(tif_list) == 0
                    
                    if force_re_extraction ~= 1
                        disp([dataset_name, 'already analysed. skipping...'])
                        continue
                    else
                    end
                else
                end
            else
            end
               
            
            disp(['Reading in avg stack for ', direc])
            
            if do_over == 1
                cd(direc)
                tif_list = dir('*.tif');
                if exist([save_path, '\tr_avg_stack.mat']) == 2
                    ave_stack = load([save_path, '\tr_avg_stack.mat']);
                    ave_stack = ave_stack.ave_stack;
                    if size(ave_stack, 3) < length(tif_list)
                        tif_start_n = size(ave_stack, 3) + 1;
                    else
                        
                        continue
                    end
                else
                    tif_start_n = 1;
                end
                
                prev_dir = pwd;
                cd(direc);
                dir_contents = dir('*.tif');
                cd(prev_dir);
                
                for tif_n = tif_start_n:length(tif_list)
                    try
                        stack_obj = ScanImageTiffReader([direc, tif_list(tif_n).name]);
                    catch
                        tif_list(tif_n) = [];
                        continue
                    end
                    curr_stack = stack_obj.data();
                    curr_stack = permute(curr_stack,[2 1 3]);
                    
                    %obtaining, logging timestamp
                    curr_time_stamp = parse_tiff_timestamp(stack_obj);
                    time_stamps(tif_n).tstamp = curr_time_stamp;
                    time_stamps(tif_n).name = dir_contents(tif_n).name;
                    
                    %checking how many color channels were acquired and saving red chan separately
                    [frame_time, zoom, n_chans, PMT_offsets, lpower, n_planes] = SI_tif_info(stack_obj);
                    
                    if extract_both_channels == 1 && n_chans == 2
                        warning('2 channels detected, extracting both channels..');
                        if del == 1
                            n_chans = 1;
                        else
                        end
                    else
                    end
                    
                    %dealing with situations where stack's frames aren't
                    %all of the same sample.
                    if n_chans == 2
                        red_stack = curr_stack(:, :, [2:2:end]);
                        curr_stack = curr_stack(:, :, [1:2:end]);
                    elseif n_planes > 1
                        curr_stack = curr_stack(:, :, [1:n_planes:end]);
                    end
                                                         
                    curr_stack = double(curr_stack);
                    ave_stack(:, :, tif_n) = std(curr_stack, 0, 3, 'omitnan');
                    disp(['Saving avg stack, tr ', int2str(tif_n), ' done.'])
                    
                    if isdir(save_path) == 0
                        mkdir(save_path);
                    else
                    end
                    
                    save([save_path, '\tr_avg_stack.mat'], 'ave_stack');
                    save([save_path, 'tif_time_stamps.mat'], 'time_stamps');
                    save([save_path, 'n_planes.mat'], 'n_planes');
                end
                
                
                
                clear ave_stack
                
            elseif do_over == 2
                
                dataset_stack = load([save_path, '\tr_avg_stack.mat']);
                dataset_stack = dataset_stack.ave_stack;
                
                ref_im = dataset_stack(:, :, 1);
                
                if ROI_type == 0
                
                    %checking if multiple ROIs have been drawn and unzipping them all
                    if exist([direc, 'ROIs\RoiSet.zip']) == 2
                        unzip([direc, 'ROIs\RoiSet.zip'], [direc, 'ROIs\']);
                    else
                    end

                    %loading in manually drawn, FIJI ROIs
                    prev_direc = pwd;
                    cd([direc, 'ROIs\']);
                    ROI_list = dir('*.roi');
                    cd(prev_direc);
                    n_ROIs = size(ROI_list, 1);
                    

                    ROI_mat = zeros(size(ref_im, 1), size(ref_im, 2), n_ROIs);
                    for ROI_n = 1:n_ROIs
                        curr_name = ROI_list(ROI_n).name;
                        curr_ROI = ReadImageJROI([direc, 'ROIs\', curr_name]);
                        ROI_mat(:, :, ROI_n) = poly2mask(curr_ROI.mnCoordinates(:, 1), curr_ROI.mnCoordinates(:, 2), size(ref_im, 1), size(ref_im, 2));
                    end
                     
                    if size(ROI_mat, 3) > 1 && fuse_ROIs == 0
                        mult_ROIs = 1;
                    else
                        mult_ROIs = 0;
                    end
                   
                elseif ROI_type == 1
                    ROI_mat = load([save_path, 'ROI_mat.mat']);
                    ROI_mat = ROI_mat.ROI_mat;
                    n_ROIs = size(ROI_mat, 3);
                    
                else
                end

                if fuse_ROIs == 1
                    ROI_mat = sum(ROI_mat, 3);
                    ROI_mat(ROI_mat > 1) = 1;
                else
                end

                %dilating ROIs if specified by user
                if dilate_ROIs > 0
                    str = strel('disk', dilate_ROIs, 0); 
                    ROI_mat = imdilate(ROI_mat, str);

                else
                end
                
                %getting rid of ROI overlap pixels if manually specified
                if remove_ROI_ovlap == 1
                    sum_mat = sum(ROI_mat, 3);
                    sum_mat = repmat(sum_mat, 1, 1, size(ROI_mat, 3));
                    ovlapi = find(sum_mat > 1);
                    ROI_mat(ovlapi) = 0;
                else
                end
                    
                
                % re-saving ROIs after fusing/dilating as specified by user
                save([save_path_base, dataset_name, '\ROI_mat.mat'], 'ROI_mat');
                
                %if this is a KC dataset, asking user to manually identify
                %ROI-warping landmarks if warping is specified
                if warp_ROIs == 1
                    
                    tiff_times = load([save_path, 'tif_time_stamps.mat']);
                    tiff_times = tiff_times.time_stamps;
                    [warping_landmarks, ROI_mat_warped] = pick_matched_landmarks(save_path);
                    
                else
                    ROI_mat_warped = [];
                    warping_landmarks = [];
                end
                
                %looping through an ave frame for each trial in the
                %dataset, for the user manually align each frame
                %to correct x-y drift, redraw ROIs to account for tissue warping, 
                %and also indicate bad z-trials
                if force_re_extr_manual == 0
                    if exist([save_path, '\ROI_mat_adj.mat']) == 2
                        continue
                    else
                    end
                else
                end
                
                done_marking = 0;
                
                try
                    n_planes = load([save_path, 'n_planes.mat']);
                    n_planes = n_planes.n_planes;
                catch
                    n_planes = 1;       %assuming n_planes is 1 for older datasets where n_planes was always 1.
                end
                
                while done_marking == 0
                    if size(ROI_mat, 3) > 1
                        if n_planes == 1
                            ROI_mat_s = sum(ROI_mat, 3);
                        elseif n_planes > 1
                            ROI_mat_s = ROI_mat(:, :, 1);
                        else
                        end
                    else
                        ROI_mat_s = ROI_mat;
                    end
%                     [lag_mat, bad_trs, done_marking, bk_ROI] = manual_xylags_zbad2(dataset_stack, ROI_mat_s, ROI_mat_warped, warping_landmarks);
                    [lag_mat, bad_trs, done_marking, bk_ROI, ROI_mat_adj] = manual_xylags_zbad2(dataset_stack, ROI_mat_s, ROI_mat_warped, [], dilate_ROIs);
                end
              
                bad_tr_list = 1:1:size(dataset_stack, 3);
                bad_tr_list(bad_trs == 1) = [];
                save([save_path, '\xy_lags.mat'], 'lag_mat');
                save([save_path, '\bad_trial_list.mat'], 'bad_tr_list');        %bad_tr_list is actually the list of good trials.
                save([save_path, '\bk_ROI.mat'], 'bk_ROI');                     %manually drawn background ROI for background subtraction.
                save([save_path, '\ROI_mat_adj.mat'], 'ROI_mat_adj');           %matrix of ROIs for each trial, including manually adjusted ROIs to account for tissue warping.
                clear lag_mat
                clear bad_trs
                
            end

        end

        
    end
end

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
        direc1 = get_path_1_direc(direc);       %in case raw direc has \1 folders created by the Suite2P program
        %checking if ROI_type is manualROIs (0) or Suite2P (1)
        if length(direc) == length(direc1)
            ROI_type = 0 ;
            save_path_base = save_path_base_manual;
        elseif length(direc) < length(direc1)
            ROI_type = 1;
            save_path_base = save_path_base_Suite2P;
            direc = direc1;
        else
        end             
        
        
        
        dataset_namei = findstr(direc, '20');
        dataset_name = direc((dataset_namei):end);
        remove_small_tifs(direc);
        prev_direc = pwd;
        cd([direc]);
        tif_list = dir('*.tif');
        
        if exist([save_path_base, dataset_name, '\trace_extraction_complete.mat']) == 2
            if force_re_extraction ~= 1
                continue
            else
            end
        else
        end
        
        disp(['Extracting traces from ' direc]);
        try
            curr_stack = ScanImageTiffReader([direc, tif_list(1).name]).data();
        catch
            keyboard
        end
        curr_stack = permute(curr_stack,[2 1 3]);
        ref_im = mean(curr_stack, 3, 'omitnan');
        
        %loading in previously saved ROI_mat
        ROI_mat = load([save_path_base, dataset_name, '\ROI_mat.mat']);
        ROI_mat = ROI_mat.ROI_mat;
        
        %if this is a KC dataset, no question of re-drawn ROIs
        if ROI_type == 0
            if size(ROI_mat, 3) > 1 && fuse_ROIs == 0
                mult_ROIs = 1;
            else
                mult_ROIs = 0;
            end
        else
        end
        
        if ROI_type == 0    %case where ROIs are manual and not Suite2P
            if mult_ROIs == 0
                ROI_mat_adj = load([save_path_base, dataset_name, '\ROI_mat_adj.mat']);
                ROI_mat_adj = ROI_mat_adj.ROI_mat_adj;
            elseif mult_ROIs == 1
                ROI_mat_adj = [];
            else
            end
        else
            ROI_mat_adj = [];
        end
        
        %extracting raw traces
        dataset_namei = findstr(direc, '20');
        dataset_name = direc((dataset_namei):end);
        
        if ROI_type == 0
            save_path_base = save_path_base_manual;
        elseif ROI_type == 1
            save_path_base = save_path_base_Suite2P;
        else
        end
        save_path = [save_path_base, dataset_name, '\' ];
        
        subtraction_spec_vec = [do_bk_subtraction, do_noisefit_subtraction];    %two manually specified booleans that specify whether or not to do the resp corrections while extracting data.
        
        %deleting old extracted traces if re-extraction is forced
        if force_re_extraction == 1
            delete([save_path, 'extracted_raw_data_mat.mat']);
        else
        end
        
        %PICK UP THREAD HERE
        %fix extract_raw_traces such that it skips every alternate frame
        %for each ROI when there are multiple planes./
        
        [raw_data_mat] = extract_raw_traces_par(direc, ROI_mat, save_path, 1, extract_both_channels, subtraction_spec_vec, force_re_extraction, ROI_mat_adj);
        
        
        %copying over files needed for further analysis to results
        %directory
        prev_direc = pwd;
        cd([direc])
        %copying over stimulus param file from raw data directory
        param_names = dir('params*.*');
        for p_file_n = 1:size(param_names, 1)
            copyfile([direc, param_names(p_file_n).name], [save_path, param_names(p_file_n).name]);            
        end
        
        % copying over PID traces to the results folder
        PID_fnames = dir('PID*.*');
        for PID_trace_n = 1:size(PID_fnames, 1)
            curr_name = PID_fnames(PID_trace_n).name;
            copyfile([direc, curr_name], [save_path, curr_name]);

        end
        %copying over first trial .tif file to results folder
        tif_fnames = dir('*.tif');
        curr_name = tif_fnames(1).name;
        copyfile([direc, curr_name], [save_path, curr_name]);
        
        disp(['Done extracting traces from ' direc]);
        cd(prev_direc)        
    end
    
end