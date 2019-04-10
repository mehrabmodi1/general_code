clear all
close all

%input variables
stimOI = 1;             %0 - only odour, 1 - odour and elec, 2 - odour and led, comma-separated vector for multiple stims
playback = 1;

%setting up initialisation variables for dataset
direc = 'C:\Data\CSHL\20150227\expt\elec_alone\';
cd(direc)
a = dir('*.tif');
n_trials = length(a);
clear a
a = dir('param*.*');
params = load(a.name);
params = params.params;
odour_list = params.odours;

if isempty(params.elec_odours) == 0
    elec_odours = str2num(params.elec_odours);
else
    elec_odours = [];
end

if isempty(params.led_odours) == 0
    led_odours = str2num(params.led_odours);
else
    led_odours = [];
end

elec_odours = [1, 2];

%building list of trials where elec stim was delivered
elec_trials = zeros(1, n_trials);
for el_odour_n = 1:length(elec_odours)
    el_odour = elec_odours(el_odour_n);
    temp = find(odour_list == el_odour);
    elec_trials(temp) = 1;    
end

elec_trials = zeros(1, n_trials) + 1;


%building list of trials where led stim was delivered
led_trials = zeros(1, n_trials);
for el_odour_n = 1:length(led_odours)
    el_odour = led_odours(el_odour_n);
    temp = find(odour_list == el_odour);
    led_trials(temp) = 1;    
end



mov_obj = VideoReader([direc 'behav_video.mj2']);           %creating movie object
n_frames = mov_obj.NumberOfFrames;      
frame1 = read(mov_obj, 1);
n_pix = size(frame1, 1).*size(frame1, 2);

figure(1)
imshow(frame1)
[x, y] = ginput(1);                     %asking user to mark the location of the fly's head
disp('Mark position of centre of fly head.');

%creating circular ROI around head 
BW = zeros(size(frame1, 1), size(frame1, 2));
BW = draw_circle(round(x), round(y), 50, BW);
BW = abs(BW - 1);                       %inverting ROI


%calculating corrcoef of each frame with frame immediately preceding it
% corr_vec = zeros(1, (n_frames - 1) );

% for frame_n = 1:(n_frames - 1)
%     if frame_n == 1
%         frame1 = single(read(mov_obj, frame_n)).*BW;
%     elseif frame_n > 1
%         frame1 = frame2;
%     end
%     
%     frame2 = single(read(mov_obj, (frame_n + 1) )).*BW;
%         
%     
%     vec1 = reshape(frame1, 1, []);
%     vec2 = reshape(frame2, 1, []);
%     
%     %calculating and storing corrcoef
%     c = corrcoef(vec1, vec2);
%     c = c(1, 2);
%     corr_vec(frame_n) = c;
%     disp(frame_n)
% end



% %building list of frame numbers in video when stimuli were delivered
int_vec = zeros(1, n_frames);
stim_frame_n = 8.*60;                   %frame number in each trial where stim is delivered (stim_time .* frames/s)
framesptrial = 14.*60;                  %no. of frames per trial
stim_frame_vec = 0:framesptrial:(framesptrial.*(n_trials - 1) );
stim_frame_vec = stim_frame_vec + stim_frame_n;


%loop to walk through each trial
response_mat = zeros(241, length(stim_frame_vec), 3) + nan;
tr_count = 0;
for stim_n = 1:length(stim_frame_vec)
    tr_count = tr_count + 1;
    %working out if current trial has elec, led or no stim
    if led_trials(stim_n) == 1 && elec_trials(stim_n) == 1
        trial_type = 3;
    elseif led_trials(stim_n) == 1
        trial_type = 2;
    elseif elec_trials(stim_n) == 1
        trial_type = 1;
    elseif led_trials(stim_n) == 0 && elec_trials(stim_n) == 0
        trial_type = 0;
    else
    end
    
    %skipping current trial if not stimOI
    if isempty(intersect(trial_type, stimOI)) == 1
        continue
    else
    end
    
    curr_stim_fr = stim_frame_vec(stim_n);
    curr_tr_frames = (curr_stim_fr - 120):(curr_stim_fr + 120);
    
    
    for frame_ni = 1:241
        frame_n = curr_tr_frames(frame_ni);
        frame = read(mov_obj, frame_n);
        
        if frame_ni == 1
            curr_tr_mat = zeros(size(frame ,1), size(frame, 2), 241);
            curr_tr_mat(:, :, 1) = frame;
        else
        end
        
        if frame_ni > 1
            curr_tr_mat(:, :, frame_ni) = frame;
        else
        end
        
         
%         frame = single(frame).*BW;
%         int_score = mean(mean(frame));
%         int_vec(1, frame_n) = int_score;
        
        disp(['done reading frame ' int2str(frame_ni)]);
    end
    
    
    %thresholding image (to get rid of still background) and calculating frame-by-frame corrcoef as a movement measure
    corr_vec = zeros(1, 241);
    diff_vec = zeros(1, 241);
    diff_vecb = zeros(1, 241);
    for curr_framen = 1:240
        if curr_framen == 1
            frame1 = curr_tr_mat(:, :, curr_framen);
            del = find(frame1 < 100);
            frame1(del) = 0;           
            frame1 = reshape(frame1, 1, []);
            
        elseif curr_framen > 1
            frame1 = frame2;
        end
        
        frame2 = curr_tr_mat(:, :, (curr_framen + 1) );
        del = find(frame2 < 100);
        frame2(del) = 0;
        frame2 = reshape(frame2, 1, []);
        
        %movt measure1       
        c = corrcoef(frame1, frame2);
        c = c(1, 2);
        corr_vec(curr_framen) = c;
        
        %movt measure2
        diff = sum(abs(frame2 - frame1))./n_pix;
        diff_vec(curr_framen) = diff;
        
        %movt measure3
        del = find(frame1 > 0);
        frame1b = frame1;
        frame1b(del) = 1;
        
        del = find(frame2 > 0);
        frame2b = frame2;
        frame2b(del) = 1;
        diffb = sum(abs(frame2b - frame1b))./n_pix;
        diff_vecb(curr_framen) = diffb; 
    end
    
    figure(2)
    subplot(2, 1, 1)    
    plot([corr_vec; diff_vecb]', '.')
    axis([1, 241, 0, 1])
    subplot(2, 1, 2)
    plot(diff_vec, '.r')
    
    response_mat(:, tr_count, 1) = corr_vec;
    response_mat(:, tr_count, 2) = diff_vec;
    response_mat(:, tr_count, 3) = diff_vecb;
    
    
    %playing video for curr trial to user
    again = 1;
    if playback == 1
        save_vid = 0;
        while again == 1
            

            for framen = 1:241
                frame = curr_tr_mat(:, :, framen);
                if save_vid == 1
                    if frame_n == 1
                        frame_mat = zeros(size(frame, 1), size(frame, 2), 241);
                        frame_mat(:, :, framen) = frame;
                    else
                        frame_mat(:, :, framen) = frame;
                    end
                      
                    
                else
                end
                
                
                
                
                figure(1)
                colormap('gray')
                imagesc(frame)
                if framen > 110 && framen < 130
                    title(['Stim ' int2str(framen) ])
                else
                    title(int2str(framen))
                end
                pause(0.01);
                drawnow
                disp(['trial no. ' int2str(stim_n)])
            end
            
            if save_vid == 1
                vid_save_dir = [direc 'saved_vid_frames\'];
                mkdir(vid_save_dir);
                for frame_n = 1:241
                    frame = frame_mat(:, :, frame_n);
                    frame = uint8(frame);
                    colormap('gray');
                    if frame_n < 150 && frame_n > 120
                        frame(10:60, 10:60) = 255;
                    else
                    end
                    imwrite(frame, [vid_save_dir 'frame' int2str(frame_n) '.jpg'], 'jpeg')
                    
                end
            else
            end
            
            again = input('Play again? 1 = yes, 0 = no')
        
            %saving video if asked to
            save_vid = input('Save Video? - 1 - yes, 0 - no')
            if save_vid == 1
                again = 1;
            else
            end
        end
        
    else
    end
    
    keyboard    
 end