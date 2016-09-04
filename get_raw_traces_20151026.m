function [raw_data_mat] = get_raw_traces_20151026(dataset, direc, saved_lags)

%extracting raw intensity data, averaged over ROI pixels and saving
%to file for further analysis
cell_roi_mat = dataset(1).ROI(2).roi;

%identifying centroids of each ROI
centroids = regionprops(cell_roi_mat, 'Centroid');
centroids = round(cat(1, centroids.Centroid));

n_cells = max(max(cell_roi_mat));
n_trials = size(dataset, 2);
curr_stack = dataset(1).imageStack;
n_frames = size(curr_stack, 3);
raw_data_mat = zeros(n_frames, n_cells, n_trials) + nan;
for trial_n = 1:n_trials
    curr_stack = dataset(trial_n).imageStack;
        
    for cell_n = 1:n_cells
        curr_stack_al = zeros(size(curr_stack, 1), size(curr_stack, 2), size(curr_stack, 3)) + nan;
        curr_lags = squeeze(saved_lags(cell_n, trial_n, :));          %This is the x,y lag by which to shift frames for this cell for this trial, as established earlier by bad_cell_catcher
        curr_lags = curr_lags';
        
        %in case current trial is a dropped trial, extracted data will be
        %dumped later on anyway, so 0 lags assumed.
        if sum(isnan(curr_lags)) > 0
            curr_lags = [0, 0];
        else
        end
        curr_offset = curr_lags;
        curr_offset(1) = curr_offset(1) .* (- 1);       %fixing sign convention
        

        %shifting all frames of current trial to extract data for this cell
        %as per lags previously identified by bad_cell_catcher
        
        if curr_offset(1) == 0 && curr_offset(2) == 0
            curr_stack_al = curr_stack;

        else
            %inserting offset corrected curr frame into curr_stack_al 
            if sign(curr_offset(1)) >= 0 && sign(curr_offset(2)) >= 0
                curr_stack_al( (1):(size(curr_stack, 1) - curr_offset(1) ), (1):(size(curr_stack, 2) - curr_offset(2) ), 1:n_frames)...
                    = curr_stack( (1+ curr_offset(1) ):(size(curr_stack, 1) ), (1+ curr_offset(2) ):(size(curr_stack, 2) ), 1:n_frames );    

            elseif sign(curr_offset(1)) < 0 && sign(curr_offset(2)) >= 0
                curr_stack_al(1 + abs(curr_offset(1) ):(size(curr_stack, 1) ), (1):(size(curr_stack, 2) - curr_offset(2) ), 1:n_frames)...
                    = curr_stack( (1):(size(curr_stack, 1) - abs(curr_offset(1)) ), (1 + curr_offset(2) ):(size(curr_stack, 2) ), 1:n_frames);

            elseif sign(curr_offset(1)) >= 0 && sign(curr_offset(2)) < 0
                curr_stack_al( (1):(size(curr_stack, 1) - curr_offset(1) ), (1 + abs(curr_offset(2)) ):(size(curr_stack, 2) ), 1:n_frames)...
                    = curr_stack(1 + (curr_offset(1) ):(size(curr_stack, 1) ), (1):(size(curr_stack, 2) - abs(curr_offset(2) ) ), 1:n_frames);

            elseif sign(curr_offset(1)) < 0 && sign(curr_offset(2)) < 0
                curr_stack_al(1 + abs(curr_offset(1) ):(size(curr_stack, 1) ), (1 + abs(curr_offset(2) ) ):(size(curr_stack, 2) ), 1:n_frames)...
                    = curr_stack( (1):(size(curr_stack, 1) - abs(curr_offset(1) ) ), (1):(size(curr_stack, 2) - abs(curr_offset(2) ) ), 1:n_frames);
            end
            
        end


        curr_roi = zeros(size(cell_roi_mat, 1), size(cell_roi_mat, 2) );
        curr_centroid = centroids(cell_n, :);
        del = find(cell_roi_mat == cell_n);
        curr_roi(del) = 1;                             %ROI of current cell
        
        for frame_n = 1:n_frames
            raw_data_mat(frame_n, cell_n, trial_n) = nanmean(nanmean(curr_stack_al(:, :, frame_n).*curr_roi) );
            
        end
        
                
    end
end