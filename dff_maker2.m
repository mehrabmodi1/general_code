function [dff_data_mat dffdata_peaks dff_bk_data_mat] = dff_maker(raw_data_mat, frame_time, Ca_width, raw_bk_data_mat)
%dff_maker syntax: dff_data_mat = dff_maker(raw_data_mat, Ca_time, frame_time)
%dff_maker takes a matrix of raw, fluorescence intensity data with frames
%along dim 1, cells along dim 2 and trials along dim 3 and returns an
%identically sized matrix with the raw intensity values normalised after
%subtracting away the baseline (dff_data_mat). dff_data_peaks is the same
%as dff_data_mat, except with peaks of width < Ca_time removed and replaced
%with baseline values. The Ca_time and frame_time inputs are
%optional, with default values of 400 and 80 ms. dff_maker also detects and
%removes peaks shorter than Ca_time determines.


%checking for Ca_time and frame_time inputs, else assigning default values
if nargin == 1
    Ca_width = 250;      %ms
    frame_time = 80;    %ms
    raw_bk_data_mat = raw_data_mat;
elseif nargin == 2
    Ca_width = 250;    %ms
    raw_bk_data_mat = raw_data_mat;
elseif nargin == 3
    raw_bk_data_mat = raw_data_mat;
end

%checking if raw_data_mat is the right size
if size(raw_data_mat, 1) == 1 | size(raw_data_mat, 2) == 1 | size(raw_data_mat, 3) == 1
    error('raw_data_mat should be a 3-D matrix with frames along dim 1, cells along dim 2 and trials along dim 3.')
else
end

%parsing raw_data_mat for information about dataset
no_frames = size(raw_data_mat, 1);
no_cells = size(raw_data_mat, 2);
no_trials = size(raw_data_mat, 3);
no_bk_frames = size(raw_bk_data_mat, 1);

%calculating dF/F using median pixel intensity as baseline
    for cell_no = 1:no_cells
        for trial_no = 1:no_trials
            dff_baseline = nanmedian(squeeze(raw_data_mat(4:no_frames, cell_no, trial_no)));
            dff_baseline_bk = nanmedian(squeeze(raw_bk_data_mat(4:no_bk_frames, cell_no, trial_no)));
            dff_data_mat(:, cell_no, trial_no) = (raw_data_mat(:, cell_no, trial_no) - dff_baseline)./dff_baseline;
            dff_bk_data_mat(:, cell_no, trial_no) = (raw_bk_data_mat(:, cell_no, trial_no) - dff_baseline_bk)./dff_baseline_bk;
        end
    end

    clear dff_baseline
    clear dff_baseline_bk
    
    
    SDs = zeros(no_cells, no_trials);
    dffdata_peaks = dff_data_mat;
    %detecting peaks in dF/F traces
    for trial_no = 1:no_trials
        for cell_no = 1:no_cells
            data_vec = squeeze(dff_data_mat(:, cell_no, trial_no));
            data_vec(1:3, 1) = 0;
            temp_SD = std(data_vec);
            %First Round of peak detection to estimate true noise SD
            if temp_SD == 0
               temp_SD = .001;
            else
            end
            [peaksi peaks] = peakdet(data_vec, temp_SD.*1.5);
            data_veci = data_vec;
            for pk_no = 1 : size(peaksi, 1)
                data_veci(uint32(peaksi(pk_no, 1)): peaksi(pk_no, 1) + 4) = NaN;
            end
            true_SD = nanstd(data_veci);
            SDs(cell_no, trial_no) = true_SD;                    %matrix of the true SDs for each cell for each trial
            clear data_veci

            %Eliminating all peaks narrower than min_width frames
            min_width = round(Ca_width./(frame_time));     %min_width is no of frames closest to being 400 ms in total    

            %thresholds for detecting initial part of Ca response and its continuation
            thresh_onset = 2.*true_SD;
            thresh = 1.*true_SD;

            data_veci = data_vec;
            pks = find(data_vec(1:(length(data_vec) - min_width))>thresh_onset);
            pks = pks';
            count = 0;
            for pk_no = 1:length(pks)
                post_pk = data_vec( (pks(1, pk_no)+1):(pks(1, pk_no)+ (min_width-1)), 1 );
                noise_points = find(post_pk < thresh);
                if isempty(noise_points) == 0
                    noise_points = noise_points  + pks(1, pk_no);
                    data_veci (noise_points) = 0;
                    data_veci (pks(1, pk_no)) = 0;
                else
                end
            end
            dffdata_peaks(:, cell_no, trial_no) = data_veci;
    
        end

    end

    clear temp_SD
    clear data_veci
    clear data_vec
    clear peaksi
    clear peaks
    clear pk_no
    clear true_SD
    clear noise_points
    