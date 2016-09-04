clear all
close all

list_path = 'D:\Data\CSHL\dataset_list_20141031.txt';
save_direc = 'D:\Data\CSHL\Analysed_data\';                             %main directory to save processed data files

norm_ratio_vec = [];
fid=fopen(list_path);
counter = 0;
while 1
    counter = counter + 1;
    direc = fgetl(fid);
    if ~ischar(direc),   break,   end
    direc = [direc '\'];
    
    
    disp(['dataset number ' int2str(counter)])
    
    
    if rem(counter, 2) == 1
        direc_w = direc;
    elseif rem(counter, 2) == 0
        direc_ag = direc;
        %----------------------------------------
        
        %reading in data for water dataset
        [n_frames, n_trials, frame_rate, zoom, d_path_w, f_name_w, tr_tg_num_w, ave_frame_w] = tiff_info(direc_w);
        [isi, odor_t_list_w, odorNames, n_reps_w, stim_time] = load_stim_param(direc_w);
        save_path_w = [save_direc d_path_w];
        frame_time = 1./frame_rate;
        stim_frame = floor(stim_time./frame_time);
        
        
        temp = sort(unique(odor_t_list_w));
        w_odor_n = temp(1, 1);
        raw_data_mati = load([save_path_w 'raw_data_mat.mat']);
        raw_data_mat = raw_data_mati.raw_data_mat;
        raw_data_mat = raw_data_mat{1, 1};

        n_cells_w = size(raw_data_mat, 2);
        clear raw_data_mati
        
        %calculating baselines for each cell, followed by dF/F for each
        %fluorescence intensity value
        Fo_mat = mean(raw_data_mat(1:(stim_frame - 1), :, :), 1);
        Fo_mat = repmat(Fo_mat, [n_frames, 1, 1]);

        dff_data_mat_w = (raw_data_mat - Fo_mat)./Fo_mat;         %calculating dF/F
        common_cell_list = load([save_path_w 'matched_cells_list.txt']);
    
    
    
        %reading in data for ag dataset
        [n_frames, n_trials, frame_rate, zoom, d_path_ag, f_name_ag, tr_tg_num_ag, ave_frame_ag] = tiff_info(direc_ag);
        [isi, odor_t_list_ag, odorNames, n_reps_ag, stim_time] = load_stim_param(direc_ag);
        save_path_ag = [save_direc d_path_ag];
        
        raw_data_mati = load([save_path_ag 'raw_data_mat.mat']);
        raw_data_mat = raw_data_mati.raw_data_mat;
        raw_data_mat = raw_data_mat{1, 1};

        n_cells_ag = size(raw_data_mat, 2);
        clear raw_data_mati
        
        %calculating baselines for each cell, followed by dF/F for each
        %fluorescence intensity value
        Fo_mat = mean(raw_data_mat(1:(stim_frame - 1), :, :), 1);
        Fo_mat = repmat(Fo_mat, [n_frames, 1, 1]);

        dff_data_mat_ag = (raw_data_mat - Fo_mat)./Fo_mat;         %calculating dF/F
        
        %--------------------------------------------
        %plotting data
        
        if common_cell_list(1, 1) == 0
            continue
        else
            n_pairs = length(common_cell_list);
            ag_odor_nx = 1:3;
            temp = find(ag_odor_nx == w_odor_n);
            ag_odor_nx(temp) = [];

            for pair_n = 1:n_pairs
                cell_pair = common_cell_list(pair_n, :);
                
                y_ax = 3;
                
                h = figure(1);
                
                plot_traces(dff_data_mat_w, odor_t_list_w, stim_frame, cell_pair(1, 1), w_odor_n, h, [2, 2, 1]);
                title(['Ringers - cell number - ' int2str(cell_pair(1, 1)) ' odor number - ' int2str(w_odor_n)]);
                axis([0, 80, -0.5, y_ax]);
                
                plot_traces(dff_data_mat_w, odor_t_list_w, stim_frame, cell_pair(1, 1), 4, h, [2, 2, 3]);
                title(['Ringers - cell number - ' int2str(cell_pair(1, 1)) ' odor number - 4' ]);
                axis([0, 80, -0.5, y_ax]);
                
                plot_traces(dff_data_mat_ag, odor_t_list_ag, stim_frame, cell_pair(1, 2), w_odor_n, h, [2, 2, 2]);
                title(['Agarose - cell number - ' int2str(cell_pair(1, 2)) ' odor number - ' int2str(w_odor_n)]);
                axis([0, 80, -0.5, y_ax]);
                
                plot_traces(dff_data_mat_ag, odor_t_list_ag, stim_frame, cell_pair(1, 2), 4, h, [2, 2, 4] );
                title(['Agarose - cell number - ' int2str(cell_pair(1, 2)) ' odor number - 4' ]);
                axis([0, 80, -0.5, y_ax]);
                
                set(gcf, 'Color', 'w');
                
                h = figure(2);
                
                plot_traces(dff_data_mat_ag, odor_t_list_ag, stim_frame, cell_pair(1, 2), ag_odor_nx(1), h, [2, 1, 1]);
                title(['Agarose - cell number - ' int2str(cell_pair(1, 2)) 'odor number - ' int2str(ag_odor_nx(1))]);
                axis([0, 80, -0.5, y_ax]);
                
                plot_traces(dff_data_mat_ag, odor_t_list_ag, stim_frame, cell_pair(1, 2), ag_odor_nx(2), h, [2, 1, 2]);
                title(['Agarose - cell number - ' int2str(cell_pair(1, 2)) 'odor number - ' int2str(ag_odor_nx(2))]);
                axis([0, 80, -0.5, y_ax]);
                
                set(gcf, 'Color', 'w');
                
                h = figure(3);
                [del, max_val_w] = plot_traces(dff_data_mat_w, odor_t_list_w, stim_frame, cell_pair(1, 1), w_odor_n, h, [2, 1, 1], 1);
                title(['Ringers - cell number - ' int2str(cell_pair(1, 1)) ' odor number - ' int2str(w_odor_n)]);
                axis([0, 80, -0.5, y_ax]);
                
                [del, max_val_ag] = plot_traces(dff_data_mat_ag, odor_t_list_ag, stim_frame, cell_pair(1, 2), w_odor_n, h, [2, 1, 2], 1);
                title(['Agarose - cell number - ' int2str(cell_pair(1, 2)) ' odor number - ' int2str(w_odor_n) ', norm ratio = ' int2str(max_val_w./max_val_ag)]);
                axis([0, 80, -0.5, y_ax]);
                set(gcf, 'Color', 'w');
                
                norm_ratio_vec = [norm_ratio_vec; max_val_w./max_val_ag];
                
                keyboard
                close all
                
            end
        
        
    
        end
    
    end
    
end
fclose(fid); 