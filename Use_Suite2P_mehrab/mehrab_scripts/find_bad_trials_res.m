function [good_trials, bad_trials] = find_bad_trials_res(direc)
%This function shortlists bad trials based on their correlation coefficients with all the other trials (post motion correction)
%and then asks the user to manually choose which of these shortlisted trials to get rid of. Trials are bad due to z-motion that could not be
%accommodated fast enough by the online z-correction algorithm.

%direc = 'D:\Data\Suite2P_registered\Suite2Ptest\20160517\1\Plane1\';
    
    %reading in trials and creating a matrix of averaged frames for each trial
    dir_contents = dir([direc, '*.tif']);
    n_trials = size(dir_contents, 1);
        
    for trial_n = 1:n_trials
        f_name = dir_contents(trial_n).name;
        
        frames = ScanImageTiffReader([direc f_name]).data();
        
        ave_frame = nanmean(frames, 3);         %calculating one averaged frame for each trial
        clear frames
        disp(['done reading in trial ' num2str(trial_n)])
        if trial_n == 1
            ave_f_mat = zeros(size(ave_frame, 1), size(ave_frame, 2), n_trials) + nan;
        else
        end
        ave_f_mat(:, :, trial_n) = ave_frame;
        
    end
    
    reshaped_mat = reshape(ave_f_mat, (size(ave_frame, 1).*size(ave_frame, 2)), n_trials);
    im_stack = ave_f_mat;
    im = reshaped_mat;
    c=corrcoef(im);                 %corrcoef matrix of mean baseline frame of each trial, across trials
    del = eye(size(c, 1) );         
    del = find(del == 1);           %getting rid of diagonal 1's
    c(del) = nan;
    clear del

    c = nanmean(c, 1);                 %vector of corrcoefs of each trial with all other trials
    corr_cutoff = .88;               %cutoff to get rid of non-matching trials just to draw ROIs with
    
    if max(c) < corr_cutoff
        corr_cutoff = max(c);
    else
    end
    
    bad_trials = find(c < corr_cutoff);
    good_trials = 1:n_trials;
    good_trials(bad_trials) = [];   %list of trials with mean corrcoef higher than cutoff defined above
    good_trials_orig = good_trials;
    bad_trials_orig = bad_trials;
    
    repeat = 1;
    repeat_counter = 1;
    while repeat == 1
        if repeat_counter == 1
            good_trials = good_trials_orig;
            bad_trials = bad_trials_orig;
        elseif repeat_counter > 1
            good_trials = 1;
            bad_trials = 2:1:n_trials;
            
        else
        end
        good_stack = mean(im_stack(:, :, good_trials), 3);
        for b_tr_n = 1:length(bad_trials)
            b_tr = bad_trials(b_tr_n);
            figure(1)
            subplot(2, 2, 1)
            imagesc(good_stack)
            colormap('gray');
            title('Mean of good trials');
            subplot(2, 2, 2)
            imagesc(im_stack(:, :, b_tr))
            title(['Bad trial # ' int2str(b_tr_n) ' of ' int2str(length(bad_trials)) ]);
            subplot(2, 2, 3)
            imagesc( abs(im_stack(:, :, b_tr) - good_stack), [min(min(good_stack)), max(max(good_stack))])
            title('Difference image');
            a = input('Keep trial? 1 - yes, 0 - no., 2 - discard dataset');
            %putting trial back into good tr list based on manual judgement call
            if a == 1
                bad_trials(b_tr_n) = 0;
                good_trials = [good_trials, b_tr];
            elseif a == 2
                b = input('Discard dataset - are you sure? 2 - Discard, 0 - Continue checking trials');
                
                if b == 2
                    del = [];
                    save([direc, 'skip_direc.txt'], 'del');
                    good_trials = [];
                    break
                else
                end
                 
                
            end

        end
        del = find(bad_trials == 0);
        bad_trials(del) = [];
        clear del
        good_tr_stack = im_stack(:, :, good_trials);
        
        repeat = input('look at badtrials again? 1 - Yes, 0 - No');
        if repeat == 1
            repeat_counter = repeat_counter + 1;
        else
        end
    end

    good_trials = sort(good_trials);
end