clear all
close all

dataset_list_paths = [%{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma.xls'};...
                      %{'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_c739_AlphaBeta.xls'}; ...
                      {'C:\Data\Data\Analysed_data\dataset_lists\dataset_list_Yoshi_PaBaEl_d5HT1b_Gamma_low_conc.xls'};...
                      ];

for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
   
    %loop to go through all experiment datasets listed in list file
    for dir_n = 1:n_dirs
        curr_dir = [dir_list{dir_n, 1}, '\'];
        cd(curr_dir);
        tiff_name = dir('*.tif');
        
        stack_obj = ScanImageTiffReader([curr_dir, tiff_name.name]);
        stack = stack_obj.data();
        stack = permute(stack,[2 1 3]);
        
        %displaying averaged frame
        figure(1)
        mean_fr = mean(stack, 3, 'omitnan');    %averaged frame for first trial
        imagesc(mean_fr, [0, 0.5.*max(max(mean_fr))]);
        %plotting on second display if present
        plot_big_fig(1);
        set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        set(gca,'ytick',[])
        set(gca,'yticklabel',[])
        colormap('gray')
        
        %loop to manually click on each visible cell
        more_cells = 1;
        ROI_mat = zeros(size(mean_fr, 1), size(mean_fr, 2));
        cell_counter = 0;
        
        %checking if cells have already been counted for this dataset
        if exist([curr_dir, 'counted_cells.mat']) == 2
            ROI_mat = load([curr_dir, 'counted_cells.mat']);
            ROI_mat = ROI_mat.ROI_mat;
            cell_counter = size(ROI_mat, 3);
        else
        end
        
        while more_cells == 1
            if cell_counter > 0
                if cell_counter > 1
                    ROIs_collapsed = max(ROI_mat, [], 3);
                else
                    ROIs_collapsed = ROI_mat;
                end
                ROIs_collapsed_rgb = ROIs_collapsed;
                ROIs_collapsed_rgb(:, :, 2:3) = zeros(size(mean_fr, 1), size(mean_fr, 2), 2);

                figure(1)
                imagesc(mean_fr, [0, 0.5.*max(max(mean_fr))]);
                set(gca,'xtick',[])
                set(gca,'xticklabel',[])
                set(gca,'ytick',[])
                set(gca,'yticklabel',[])
                
                
                %adding on a layer of markers for cells already counted.
                hold on
                h = imagesc(ROIs_collapsed_rgb);
                h.AlphaData = ROIs_collapsed;
                hold off
                
            else
            end
            
            title(['Cell count: ' num2str(cell_counter)])
            
            [x,y] = ginput(1);
            
            %checking of click was out of bounds to initialise end-count question
            if x < 0 || x > size(mean_fr, 2) || y < 0 || y > size(mean_fr, 1)
                choice = questdlg('Done counting cells?', 'Cell-counting', 'Done', 'Continue', 'Remove previous cell', 'Continue');
                if strcmp(choice, 'Done') == 1
                    more_cells = 0;
                elseif strcmp(choice, 'Continue') == 1
                    more_cells = 1;
                elseif strcmp(choice, 'Remove previous cell') == 1
                    if cell_counter > 0
                        ROI_mat(:, :, cell_counter) = [];
                        cell_counter = cell_counter - 1;
                    else
                    end
                end
                
            else
                cell_counter = cell_counter + 1;
                curr_ROI = zeros(size(mean_fr, 1), size(mean_fr, 2));
                curr_ROI = draw_circle(y, x, 4, curr_ROI, 1);
                ROI_mat(:, :, cell_counter) = curr_ROI;
            end
            
            %saving counted ROI matrix every 5 cells
            if rem(cell_counter, 5) == 0
                save([curr_dir, 'counted_cells.mat'], 'ROI_mat');
            else
            end
        end
        save([curr_dir, 'counted_cells.mat'], 'ROI_mat');
        remote_path = '\\dm11\turnerlab\Mehrab\Analysed_data';
        curr_dir_namei = findstr(curr_dir, '\20');
        curr_dir_name = curr_dir(curr_dir_namei:end);
        save([remote_path, curr_dir_name, 'counted_cells.mat'], 'ROI_mat')
               
    end
    
end
close figure 1