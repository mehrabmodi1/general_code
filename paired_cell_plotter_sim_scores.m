clear all
close all

list_path = 'C:\Data\CSHL\dataset_list_20141031.txt';
save_direc = 'C:\Data\CSHL\Analysed_data\';                             %main directory to save processed data files

norm_ratio_vec = [];
saved_corrs = [];
saved_rmses = [];
saved_corrs_r = [];
saved_rmses_r = [];
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
                
                
                %-------------------------------------
                %calculating similarity between pair ave dF/F curves as
                %well as randomly chosen cells as controls.
                odor_w = unique(odor_t_list_w);                                         %odor number used in the ringer's dataset
                odor_w = odor_w(1, 1);
                t_list_w = find(odor_t_list_w == odor_w);
                t_list_ag = find(odor_t_list_ag == odor_w);
                
                cell_traces_w = dff_data_mat_w( (stim_frame - 10):end, cell_pair(1, 1), t_list_w);
                cell_traces_ag = dff_data_mat_ag( (stim_frame - 10):end, cell_pair(1, 2), t_list_ag);
                mean_trace_w = nanmean(cell_traces_w, 3);
                mean_trace_w = mean_trace_w./max(mean_trace_w);
                mean_trace_ag = nanmean(cell_traces_ag, 3);
                mean_trace_ag = mean_trace_ag./max(mean_trace_ag);
                
                c = corrcoef(mean_trace_w, mean_trace_ag);
                c = c(1, 2);
                saved_corrs = [saved_corrs; c];
                
                rmse = sqrt(mean((mean_trace_w - mean_trace_ag).^2));
                saved_rmses = [saved_rmses; rmse];
                
                
                %calculating similarity scores for water and ag traces with
                %those of 30 randomly picked cells from the other dataset
                rand_cell_list_w = round( (rand(30, 1).*(n_cells_ag - 1)) ) + 1;
                c_vec_r = zeros(30, 1) + nan;
                rmse_vec_r = zeros(30, 1) + nan;
                for r_cell_n = 1:30
                    r_cell = rand_cell_list_w(r_cell_n);
                    r_traces = dff_data_mat_ag( (stim_frame - 10):end, r_cell, t_list_ag);
                    r_trace_ave = nanmean(r_traces, 3);
                    r_trace_ave = r_trace_ave./max(r_trace_ave);
                    
                    c_r = corrcoef(mean_trace_w, r_trace_ave);
                    c_vec_r(r_cell_n) = c_r(1, 2);
                    
                    rmse_r = sqrt(mean((mean_trace_w - r_trace_ave).^2));
                    rmse_vec_r(r_cell_n) = rmse_r;
                   
                end
                saved_corrs_r = [saved_corrs_r; mean(c_vec_r)];
                saved_rmses_r = [saved_rmses_r; mean(rmse_vec_r)];
                
                
                rand_cell_list_ag = round( (rand(30, 1).*(n_cells_w - 1)) ) + 1;
                c_vec_r = zeros(30, 1) + nan;
                rmse_vec_r = zeros(30, 1) + nan;
                for r_cell_n = 1:30
                    r_cell = rand_cell_list_ag(r_cell_n);
                    r_traces = dff_data_mat_w( (stim_frame - 10):end, r_cell, t_list_w);
                    r_trace_ave = nanmean(r_traces, 3);
                    r_trace_ave = r_trace_ave./max(r_trace_ave);
                    
                    c_r = corrcoef(mean_trace_w, r_trace_ave);
                    c_vec_r(r_cell_n) = c_r(1, 2);
                    
                    rmse_r = sqrt(mean((mean_trace_w - r_trace_ave).^2));
                    rmse_vec_r(r_cell_n) = rmse_r;
                   
                end
                saved_corrs_r = [saved_corrs_r; mean(c_vec_r)];
                saved_rmses_r = [saved_rmses_r; mean(rmse_vec_r)];
                
                
                
                
            end
        
        
    
        end
    
    end
    
end
fclose(fid); 

%plotting ave curve similarity scores
figure(1)
means_corr = [mean(saved_corrs), mean(saved_corrs_r)];
ses_corr = [std(saved_corrs), std(saved_corrs_r)];
n_vec = sqrt([length(saved_corrs), length(saved_corrs_r)]);
ses_corr = ses_corr./n_vec;
barweb(means_corr', ses_corr')
ylabel('Correlation coefficients')
hold on

r_vec1 = rand(1, length(saved_corrs)).*0.25 + .75;
plot(r_vec1, saved_corrs, 's', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'b')
hold on
r_vec2 = rand(1, length(saved_corrs_r)).*0.25 + 1.01;
plot(r_vec2, saved_corrs_r, 's', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'r')
set(gcf, 'Color', 'w')


figure(2)
means_corr = [mean(saved_rmses), mean(saved_rmses_r)];
ses_corr = [std(saved_rmses), std(saved_rmses_r)];
n_vec = sqrt([length(saved_rmses), length(saved_rmses_r)]);
ses_corr = ses_corr./n_vec;
barweb(means_corr', ses_corr')
ylabel('Root Mean Square Error')
hold on

r_vec1 = rand(1, length(saved_rmses)).*0.25 + .75;
plot(r_vec1, saved_rmses, 's', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'b')
hold on
r_vec2 = rand(1, length(saved_rmses_r)).*0.25 + 1.01;
plot(r_vec2, saved_rmses_r, 's', 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'r')
set(gcf, 'Color', 'w')

[h1, p1] = ttest2(saved_corrs, saved_corrs_r);
[h2, p2] = ttest2(saved_rmses, saved_rmses_r);
