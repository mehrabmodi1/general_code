function [dff_data_mat dff_data_mat_orig dffdata_peaks dffdata_peaks_orig trial_type_vec bad_cellsi bad_trialsi dff_data_mat_nans] = imaging_data_sanitiser(dff_data_mat, dffdata_peaks, trial_type_vec)
%imaging_data_sanitiser syntax: [dff_data_mat dff_data_mat_orig trial_type_vec bad_cellsi] = imaging_data_sanitiser(dff_data_mat, trial_type_vec)
%This function takes a dff_data_mat matrix and eliminates especially bad
%trials, and even slightly bad cells. An eliminated cell is eliminated from
%all trials. Scoring for 'badness' is based on the number of infs and nans
%in the respective data sub-set.


no_frames = size(dff_data_mat, 1);
no_cells = size(dff_data_mat, 2);
no_trials = size(dff_data_mat, 3);
dff_data_mat_orig = dff_data_mat;
dffdata_peaks_orig = dffdata_peaks;
dff_data_mat_nans = dff_data_mat;

    %detecting large numbers of infs/nans in dff data, over entire dataset
    infs = sum(sum(sum(isinf(dff_data_mat))));
    nans = sum(sum(sum(isnan(dff_data_mat))));
    fraction = (infs + nans)./(size(dff_data_mat, 1).*size(dff_data_mat, 2).*size(dff_data_mat, 3));
    
    %only analysing datasets with significant nos of errors further
    if  infs > 100 | nans > 100
    
       %detecting large numbers of infs/nans in dff data, on a trial by trial
       %basis
       
       t_scores = zeros(no_trials, 2);
       for trial_no = 1:no_trials
           data = squeeze(dff_data_mat(:, :, trial_no));
           infs_t = sum(sum(isinf(data)));
           nans_t = sum(sum(isnan(data)));
           t_scores(trial_no, :) = [infs_t, nans_t];

       end


       %LOOKING FOR BAD TRIALS
       %identifiying individual trials that are bad (not rescuable) as
       %those with very high error scores. these trials to be eliminated
       %entirely, keeping track of eliminated trials in trial_type vec.
       bad_trialsi1 = find(t_scores(:, 1)> 2000);
       bad_trialsi2 = find(t_scores(:, 2)> 2000);
       bad_trialsi = unique([bad_trialsi1; bad_trialsi2]);
       clear bad_trialsi1
       clear bad_trialsi2
              
       
       dff_data_mat_nans(:, :, bad_trialsi) = nan;        %eliminating bad trials
             
       
       dff_data_mat(:, :, bad_trialsi) = [];        %eliminating bad trials
       dffdata_peaks(:, :, bad_trialsi) = [];
       trial_type_vec(bad_trialsi) = [];
       
       %generating a warning message if > 5% of trials are bad
       if length(bad_trialsi) > round(no_trials./20)
           disp(['WARNING, more than 5% trials irretrievably bad (discarded)']);
           
       else
       end
       
       %Counting number of affected trails to see if rescue by elimination
       %of cells is needed. these cells to be eliminated from entire dataset.
       affected_trialsi1 = find(t_scores(:, 1) > 100);
       affected_trialsi2 = find(t_scores(:, 2) > 100);
       affected_trials = unique([affected_trialsi1; affected_trialsi2]);
       length(affected_trials);
       no_trialsx = size(dff_data_mat, 3);
       %keyboard
       %LOOKING FOR BAD CELLS
       %condition for looking for bad cells is that more than 5% of trials
       %should have some  infs/nans
      
       if length(affected_trials) > round(no_trials./20)
           saved_cell_scores = zeros(2, no_cells, no_trials);  

           for trial_no = 1:no_trialsx
               for cell_no = 1:no_cells
                   data = squeeze(dff_data_mat(:, cell_no, trial_no));
                   infs = sum(sum(isinf(data)));
                   nans = sum(sum(isnan(data)));
                   saved_cell_scores(:, cell_no, trial_no) = [infs, nans];

               end
           end

           
           %bad cell criterion: if a cell has more nans than 1/3 the no of
           %frames in a trial over the entire dataset; it is dropped. 
           [null, bad_cellsi] = find(saved_cell_scores(2, :, :) > round(no_frames.*.3));   
           bad_cellsi = unique(rem(bad_cellsi, no_cells));
           temp = find(bad_cellsi == 0);
           bad_cellsi(temp) = no_cells;
          
           %dropping bad cells from dff_data_mat
           dff_data_mat(:, bad_cellsi, :) = [];
           dffdata_peaks(:, bad_cellsi, :) = [];
       else
           bad_cellsi = [];
           bad_trialsi = [];
                      
       end
       
    else
        bad_cellsi = [];
        bad_trialsi = [];

    end   
   
    %checking for excessive depletion of dataset
    if (size(dff_data_mat, 1).* size(dff_data_mat, 2).* size(dff_data_mat, 1)) < (size(dff_data_mat_orig, 1).* size(dff_data_mat_orig, 2).* size(dff_data_mat_orig, 1).* 0.5)
        disp(['WARNING; more than half the dataset has been dropped. Ratio- ' num2str((size(dff_data_mat, 1).* size(dff_data_mat, 2).* size(dff_data_mat, 1))./(size(dff_data_mat_orig, 1).* size(dff_data_mat_orig, 2).* size(dff_data_mat_orig, 1)))])
                
    else 
    end
  