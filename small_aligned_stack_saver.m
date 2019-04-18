save_refim = 1;

%loading first 15 frames of last stack
direc = 'D:\Data\CSHL\Resonant\movt_imgs\';
old_dir = pwd;
cd(direc);
dir_contents = dir('*.tif');
n_files = size(dir_contents, 1);

%identifying newest stack
for file_n = 1:n_files
    curr_date = dir_contents(file_n).date;
    curr_date = datenum(curr_date);
    if file_n == 1
        max_date = curr_date;
        max_n = file_n;
        continue
    else
        if (curr_date < max_date) == 1
            continue
        elseif (curr_date > max_date) == 1
            max_date = curr_date;
            max_n = file_n;             %last file_n
        end
    end
end
last_file = dir_contents(max_n).name;

%loading lateral movt lags saved by Scanimage
f_name_text = last_file(1:(end - 9));
f_name_rest = last_file((end - 9):(end - 3));
lag_f_name = [f_name_text, 'Motion', f_name_rest, 'csv'];
[del,del1,raw] = xlsread(lag_f_name);
lags = raw(2:end, 5);           %lags stored as a cell array

n_frames = 40;    
for frame_n = 1:n_frames
    curr_frame = imread(last_file, frame_n);
    curr_lags = lags{frame_n};
    space_i = findstr(curr_lags, ' ');
    xlag = str2num(curr_lags(2:space_i));
    ylag = str2num(curr_lags((space_i + 1):(end - 1)));
    curr_frame_s = shiftmatrix(curr_frame, 1, [(xlag.*-1), (ylag.*-1)], nan);
    
    if frame_n == 1
        ave_frames = zeros(size(curr_frame, 1), size(curr_frame, 2), n_frames) + nan;
        ave_frames(:, :, 1) = double(curr_frame_s);
    else
        ave_frames(:, :, frame_n) = double(curr_frame_s); 
    end
end
ave_frame = nanmean(ave_frames, 3);

if save_refim ==1
    save([last_file(1:(end - 4)) '_ave.mat'], 'ave_frame')
else
end

cd(old_dir)

