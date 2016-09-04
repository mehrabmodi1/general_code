clear all
close all


direc_lists_mat = [{'D:\Data\CSHL\dataset_list_sustained_MB418B_20160414.txt'};... %KC A'B'
                    {'D:\Data\CSHL\dataset_list_sustained_MB185B_20160426.txt'};... %KC AB
                    %{'D:\Data\CSHL\dataset_list_stim_dur_20160316.txt'} ... KCs all
                    {'D:\Data\CSHL\dataset_list_sustained_MB131B_20160515.txt'}... %KC G
                    ]; 

                    
color_vec = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\std_color_vec.txt');
odor_name_list = load('C:\Users\Mehrab\Google Drive\Backup\Stuff\CSHL\Glenn lab\Code\odor_name_list.mat');
odor_name_list = odor_name_list.odor_name_list;

a = colormap('bone');
greymap = flipud(a);
colormap(greymap)

%odor_list = [3, 6];
odor_list= [1, 2, 3, 5, 7];


suppress_plots = 1;       %1 - doesn't plot stuff, 0 - plots stuff

%figure property initialisation variables
plot_height = 200;
plot_width = 200;
axis_font_size = 8;
tick_length = [0.02, 0.05];
line_widths = 2;
marker_sizes = 1;
marker_sizes_f = 4;
box_l_width = 0.5;


%%                    
n_direc_lists = size(direc_lists_mat, 1);
ratio_mat_saved = [];
ratio_mat_sd_saved = [];
ratio_mat_saved_col = [];
ratio_mat_saved_col_sd = [];

for direc_list_n = 1:n_direc_lists
    list_direc = direc_lists_mat{direc_list_n, 1};
    fid = fopen(list_direc);
    direc_counter = 0;
        
    %loop to go through all experiment datasets listed in list file
    while 1
        direc_counter = direc_counter + 1;
        direc = fgetl(fid);

        if ischar(direc) ~= 1
            break
        else
        end

        direc = [direc, '\'];

        %replacing C: with D:
        a = strfind(direc, 'C:');
        direc(a) = 'D';

        %skipping a particular OK107 dataset where n_frames was < stim time
        del = strfind(direc, '20160310');
        if isempty(del) == 0
            direc_counter = direc_counter - 1;
            continue
        else
        end
        

        n_short_selective = 0;
        n_long_selective = 0;
        
        try
            curr_data = load([direc 'cell_classifications.mat']);
        catch
            disp('dataset missing...')
            continue
        end
            
        curr_data = curr_data.classific_data;

        if direc_counter == 1
            length_vecs = size(curr_data.r_vecs, 1);
            n_odors = size(curr_data.r_vecs, 3);

            all_r_vecs = curr_data.r_vecs;
            all_h_vecs = curr_data.h_vecs;
        else
            n_cells_old = size(all_r_vecs, 2);
            n_cells_new = size(curr_data.r_vecs, 2);

            all_r_vecsi = all_r_vecs;
            all_r_vecs = zeros(length_vecs, (n_cells_old + n_cells_new), n_odors);
            all_r_vecs(:, 1:n_cells_old, :) = all_r_vecsi;
            all_r_vecs(:, (n_cells_old + 1):(n_cells_old + n_cells_new), :) = curr_data.r_vecs;

            all_h_vecsi = all_h_vecs;
            all_h_vecs = zeros(length_vecs, (n_cells_old + n_cells_new), n_odors);
            all_h_vecs(:, 1:n_cells_old, :) = all_h_vecsi;
            all_h_vecs(:, (n_cells_old + 1):(n_cells_old + n_cells_new), :) = curr_data.h_vecs;
        end
        
        
    end
    %checking for absence of 1s stim and appearance of on or sus win
    %responses or vice versa
    h_vecs = all_h_vecs;

    %getting rid of cells with missing data (either 1s or all longer response data)
    nan_mat = isnan(h_vecs);
    del = squeeze(nan_mat(6, :, :));
    temp = find(del == 1);
    [nan_1s_cells, del] = ind2sub([size(del, 1), size(del, 2)], temp );

    h_vecsi = h_vecs;       %changing this copy of h_vecs to only consider cells for whom both, 1s and longer stim responses are available
    h_vecsi(:, nan_1s_cells, :) = nan;

    del = squeeze(nan_mat(3, :, :));
    temp = find(del == 1);
    [nan_long_cells, del1] = ind2sub([size(del, 1), size(del, 2)], temp);
    nan_long_cells = unique(nan_long_cells);
    h_vecsi(:, nan_long_cells, :) = nan;



    n_cells = size(h_vecs, 2);
    for odor_n = 1:n_odors
        short_ON_list = find(h_vecsi(6, :, odor_n) == 1);
        short_OFF_list = 1:n_cells;
        short_OFF_list(short_ON_list) = [];

        long_ON_list = find(h_vecsi(3, :, odor_n) == 1);
        long_OFF_list = 1:n_cells;
        long_OFF_list(long_ON_list) = [];

        %1. 1s ON AND longer OFF
        short_selective_list = intersect(short_ON_list, long_OFF_list);

        %2. 1s OFF AND longer ON
        long_selective_list = intersect(short_OFF_list, long_ON_list);

        %3. 1s ON AND longer ON
        both_on_list = intersect(short_ON_list, long_ON_list);

        %4. 1s ON OR longer ON
        responder_list = union(short_ON_list, long_ON_list);

        if odor_n == 1
            n_short_on = 0;
            n_long_on = 0;
            n_non_selective_on = 0;
            n_short_selective = 0;
            n_long_selective = 0;
            n_responders = 0;
            n_cells_saved = 0;

        else

            n_short_on = n_short_on + length(short_ON_list);                            %summed across odors
            n_long_on = n_long_on + length(long_ON_list);                               %summed across odors
            n_non_selective_on = n_non_selective_on + length(both_on_list);             %summed across odors
            n_short_selective = n_short_selective + length(short_selective_list);       %summed across odors
            n_long_selective = n_long_selective + length(long_selective_list);          %summed across odors
            n_responders = n_responders + length(responder_list);
            
        end
    end
    
    saved_short_selectives(direc_counter, direc_list_n) = n_short_selective./(n_cells_saved.*n_odors);
    saved_long_selectives(direc_counter, direc_list_n) = n_long_selective./(n_cells_saved.*n_odors);
%%
    n_odors = size(all_r_vecs, 3);
    n_cells = size(all_r_vecs, 2);
    n_vec_elements = size(all_r_vecs, 1);

    %Plotting summed count distributions of each type of peak for all odours
    ave_h_vecs = nanmean(all_h_vecs, 2);
    %Looking at each cell individually
    unresp_h_vec = zeros(n_odors, 1);
    for odor_n = 1:n_odors
        unresp_counter = 0;
        for cell_n = 1:n_cells
            curr_h_vec = all_h_vecs(:, cell_n, odor_n);
            if sum(curr_h_vec) == 0
                unresp_counter = unresp_counter + 1;

            else
            end
        end
        unresp_h_vec(odor_n, 1) = unresp_counter;
    end
    ave_h_vecs(6, 1, :) = unresp_h_vec./n_cells;
    ave_h_vecs_saved(:, :, direc_list_n) = squeeze(ave_h_vecs);     %keeping track of cell type vecs to compare across dataset lists
    
    figure(1)
    for odor_n = 1:n_odors
        odor_ni = odor_list(odor_n);
        curr_color = color_vec(odor_ni, :);
        curr_odor_name = odor_name_list{odor_ni, 1};
        plot(ave_h_vecs(:, odor_n), 'LineWidth', 2, 'Color', curr_color)
        odor_name_vec{odor_n, 1} = curr_odor_name;
        hold on
        
    end
    legend(odor_name_vec, 'Location','eastoutside','Orientation','vertical');
    title('Probability of sig response of each type')
    ylabel('Pobability of detection of significant response')
    ax = gca;
    ax.XTick = [1:6];
    ax.XTickLabel = {'On', 'Ramp', 'Sus', 'Off', 'Multi-pk', 'Unresp'};


    %Plotting response size distributions of each type of peak for all odours
    all_r_vecsi = all_r_vecs;
    del = find(all_r_vecsi == 0);
    all_r_vecsi(del) = nan;
    ave_r_vecs = nanmean(all_r_vecs, 2);
    sd_r_vecs = nanstd(all_r_vecs, [], 2);

    figure(2)
    for odor_n = 1:n_odors
        odor_ni = odor_list(odor_n);
        curr_color = color_vec(odor_ni, :);
        shadedErrorBar([], ave_r_vecs(:, odor_n), sd_r_vecs(:, odor_n), {'-o', 'markerfacecolor', curr_color}, 1);
        hold on
    end
    title('Mean +/- SD response size of each type, across odors')
    ylabel('Significant response size (dF/F)')
    ax = gca;
    ax.XTick = [1:4];
    ax.XTickLabel = {'On', 'Ramp', 'Sus', 'Off'};


    %making 3D plot of response size along 3 of the windows
    figure(3)
    for odor_n = 1:n_odors
        odor_ni = odor_list(odor_n);
        curr_color = color_vec(odor_n, :);
        plot3(all_r_vecs(1, :, odor_n), all_r_vecs(3, :, odor_n), all_r_vecs(4, :, odor_n), 'LineStyle', 'none', 'Color', curr_color, 'Marker', 'O', 'MarkerFaceColor', curr_color)
        hold on
        xlabel('On response (dF/F)')
        ylabel('Sus response (dF/F)')
        zlabel('Off response (dF/F)')
    end
        
    %calculating ratio of on window response size to sus window response
    %size for all significantly responsive cells in either the on or sus
    %windows.
    ratio_mat = squeeze(all_r_vecs(1, :, :) ./ all_r_vecs(3, :, :));
    ratio_mat_saved = [ratio_mat_saved; nanmedian(ratio_mat, 1)];
    ratio_mat_sd_saved = [ratio_mat_sd_saved; mad(ratio_mat, 1, 1)];
    
    ratio_mat_saved_col = [ratio_mat_saved_col; nanmedian(reshape(ratio_mat, [], 1))];
    ratio_mat_saved_col_sd = [ratio_mat_saved_col_sd; mad(reshape(ratio_mat, [], 1), 1, 1)'];
       
    
end

%%
%Plotting ave_h_vecs across dataset lists
figure(5)
for direc_list_n = 1:n_direc_lists
    curr_h_vecs = ave_h_vecs_saved(:, :, direc_list_n);
    o_ave_h_vec = nanmean(curr_h_vecs, 2);                  %h_vec averaged across odours
    
    plot(o_ave_h_vec, 'LineWidth', 2)
    hold on
end
set_name_vec = {'ApBp', 'AB', 'G'};
legend(set_name_vec, 'Location','eastoutside','Orientation','vertical');
title('Probability of sig response of each type')
ylabel('Pobability of detection of significant response')
ax = gca;
ax.XTick = [1:6];
ax.XTickLabel = {'On', 'Ramp', 'Sus', 'Off', 'Multi-pk', 'Unresp'};



figure(6)
plot(ratio_mat_saved', 'LineWidth', 2)
hold on
legend({'ApBp', 'AB', 'G'}, 'Location','eastoutside','Orientation','vertical');
title('Ratio of On resp to Sus resp')
ylabel('On response ratio')
ax = gca;
ax.XTick = [1:5];
%ax.XTickLabel = {'On', 'Ramp', 'Sus', 'Off', 'Multi-pk', 'Unresp'};
ax.XTickLabel = odor_name_vec;


figure(7)
for direc_list_n = 1:n_direc_lists
    shadedErrorBar([], ratio_mat_saved(direc_list_n, :), ratio_mat_sd_saved(direc_list_n, :), {'LineWidth', 2}, 1)
    hold on
end
legend({'ApBp', 'AB', 'G'}, 'Location','eastoutside','Orientation','vertical');
title('Ratio of On resp to Sus resp')
ylabel('On response ratio')
ax = gca;
ax.XTick = [1:5];
%ax.XTickLabel = {'On', 'Ramp', 'Sus', 'Off', 'Multi-pk', 'Unresp'};
ax.XTickLabel = odor_name_vec;



figure(8)
bw_legend = [];
errorbar(ratio_mat_saved_col', ratio_mat_saved_col_sd')
title('Ratio of On resp to Sus resp')
ylabel('On response ratio')
ax = gca;
ax.XTick = [1:5];
ax.XTickLabel = {'ApBp', 'AB', 'G'};

