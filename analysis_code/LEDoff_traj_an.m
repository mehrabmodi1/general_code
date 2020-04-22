clear all
close all

base_path = 'C:\Data\Data\Analysed_data\data_sharing\Data_from_Yoshi\upwind_screen_trajectories\';
fname_empty = 'Results_10sLED-air_Empty-split';
fname_MB077B = 'Results_10sLED-air_MB077B';

frame_time = 1./30; %in s

%reading in parsed data
[ctrl_data_mat, ctrl_params, ctrl_fnums, MBON_data_mat, MBON_params, MBON_fnums, LED_oni, LED_offi] = load_data(base_path, fname_empty, fname_MB077B, frame_time);
%getting distances from arena center
ctrl_c_dists = dist_from_center(ctrl_data_mat, ctrl_params);
keyboard

%----------------
%worker functions
function [ctrl_data_mat, ctrl_params, ctrl_fnums, MBON_data_mat, MBON_params, MBON_fnums, LED_oni, LED_offi] = load_data(base_path, fname_empty, fname_MB077B, frame_time)

if exist([base_path, fname_empty, '_data.mat']) ~= 1
    
    %reading in and parsing trajectory data from file
    [ctrl_data_mat, ctrl_meta_data] = read_data_file([base_path, fname_empty, '.txt']);
    [MBON_data_mat, MBON_meta_data] = read_data_file([base_path, fname_MB077B, '.txt']);
    
    %writing sorted, parsed data to .mat files for quicker loading 
    save([base_path, fname_empty, '_data.mat'], 'ctrl_data_mat');
    save([base_path, fname_empty, '_metdata.mat'], 'ctrl_meta_data');
    save([base_path, fname_MB077B, '_data.mat'], 'MBON_data_mat');
    save([base_path, fname_MB077B, '_metdata.mat'], 'MBON_meta_data');
else
    %loading in previously parsed data
    ctrl_data_mat = load([base_path, fname_empty, '_data.mat']);
    ctrl_data_mat = ctrl_data_mat.ctrl_data_mat;
    ctrl_meta_data = load([base_path, fname_empty, '_metdata.mat']);
    ctrl_meta_data = ctrl_meta_data.ctrl_meta_data;
    MBON_data_mat = load([base_path, fname_MB077B, '_data.mat']);
    MBON_data_mat = MBON_data_mat.MBON_data_mat;
    MBON_meta_data = load([base_path, fname_MB077B, '_metdata.mat']);
    MBON_meta_data = MBON_meta_data.MBON_meta_data;
end
ctrl_fnums = ctrl_data_mat(1, :);
MBON_fnums = MBON_data_mat(1, :);
LED_oni = find(ctrl_fnums(1, :) == 0);      %frame_n when LED came on; checked that it's the same for ctrl and MBON datasets
LED_offi = LED_oni + round(10./frame_time);

ctrl_data_mat(1, :) = [];
MBON_data_mat(1, :) = [];
ctrl_data_mat = sortrows(ctrl_data_mat);
MBON_data_mat = sortrows(MBON_data_mat);
ctrl_params = ctrl_data_mat(:, 1:5);        %timestamp, arena_n, cam_n, fly_n and param_n for each row of data_mat
MBON_params = MBON_data_mat(:, 1:5);        %timestamp, arena_n, cam_n, fly_n and param_n for each row of data_mat
ctrl_data_mat(:, 1:5) = [];
MBON_data_mat(:, 1:5) = [];

function [data_mat, met_data_mat] = read_data_file(file_path)
%This function reads in relevant data as doubles in data mat and metadata
%in met_data_mat. The first five columns of data_mat are time stamp, arena_n, cam_n, fly
%number and parameter number (ie. X, Y or Heading) in that order
fid = fopen(file_path);
curr_line = 0;
data_mat = [];
met_data_mat = [];
line_n = 0;
while 1
    line_n = line_n + 1;
    curr_line = fgetl(fid);
    if curr_line == -1
        break
    else
    end
    if line_n == 1
        sub_str = curr_line(8);
    else
    end
    
    delimiteri = findstr(curr_line, sub_str);
    curr_row = cell(1, 10);
    %reading in metadata text
    for col_n = 1:10
        curr_delimiters = [(delimiteri(col_n) + 1), delimiteri((col_n + 1))];        
        curr_num = str2num(curr_line(curr_delimiters(1):curr_delimiters(2)));
        if isempty(curr_num) == 1
            curr_num = {curr_line(curr_delimiters(1):curr_delimiters(2))};
        else
        end
        curr_row{1, col_n} = curr_num;
    end
    met_data_mat = [met_data_mat; curr_row];    %metadata matrix
    
    %reading in numeric data
    curr_row = zeros(1, (length(delimiteri) - 10));
    for col_n = 11:length(delimiteri)
        if col_n < length(delimiteri)
            curr_delimiters = [(delimiteri(col_n) + 1), delimiteri((col_n + 1))];
        elseif col_n == length(delimiteri)
            curr_delimiters = [(delimiteri(col_n) + 1), length(curr_line)];
        else
        end
        
        curr_num = str2num(curr_line(curr_delimiters(1):curr_delimiters(2)));
        curr_row(1, (col_n - 10)) = curr_num;        %data matrix
            
    end
    %adding t_stamp, fly_n and param_type (x - 0, y - 1, heading - 2) in that order at row beginning to allow easy sorting
    if line_n > 1
        %1. tstamp
        col_n = 9;
        curr_delimiters = [(delimiteri(col_n) + 1), (delimiteri(col_n + 1) - 1)];
        tstamp = curr_line(curr_delimiters(1):curr_delimiters(2));
        tstamp(9) = [];
        tstamp = str2double(tstamp).*1e-13; %mapping time stamp (yyyymmddhhmmss) down to reasonable magnitude
        
        %2. Arena n (1 or 2)
        col_n = 7;
        curr_delimiters = [(delimiteri(col_n) + 1), (delimiteri(col_n + 1) - 1)];
        arena_n = curr_line(curr_delimiters(1):curr_delimiters(2));
        if strcmp(arena_n, 'Arena2') == 1
            arena_n = 2;
        elseif strcmp(arena_n, 'olfactoryArena1') == 1 
            arena_n = 1;
        else
            keyboard
        end
               
        %3. Cam n (0 0r 1)
        col_n = 8;
        curr_delimiters = [(delimiteri(col_n) + 1), (delimiteri(col_n + 1) - 1)];
        cam_n = curr_line(curr_delimiters(1):curr_delimiters(2));
        if strcmp(cam_n, 'Cam0') == 1
            cam_n = 0;
        elseif strcmp(cam_n, 'Cam1') == 1 
            cam_n = 1;
        else
            keyboard
        end
        
        %4. fly_n
        col_n = 10;
        curr_delimiters = [(delimiteri(col_n) + 1), (delimiteri(col_n + 1) - 1)];
        fly_n = str2num(curr_line(curr_delimiters(1):curr_delimiters(2)));
        
        %5. param_n
        col_n = 2;
        curr_delimiters = [(delimiteri(col_n) + 1), (delimiteri(col_n + 1) - 1)];
        param_n = curr_line(curr_delimiters(1):curr_delimiters(2));
        if strcmp(param_n, 'X') == 1
            param_n = 0;
        elseif strcmp(param_n, 'Y') == 1
            param_n = 1;
        elseif strcmp(param_n, 'Heading') == 1
            param_n = 2;
        else
        end
        params = [tstamp, arena_n, cam_n, fly_n, param_n];
    elseif line_n == 1
        params = zeros(1, 5) + nan;  %header row
    else
    end
    curr_row = [params, curr_row];
    data_mat = [data_mat; curr_row];
       
end
fclose(fid);
end
end

function [param_mat, data_mat] = get_good_flies(param_mat, data_mat)
%This function hierarchically searches param mat for three rows of
%measurements for the same fly, for each fly. It then concatenates such
%good row triplets together to generate a cleaned dataset of observations
%for further analysis.

t_list = unique(param_mat(:, 1));   %list of unique time stamps
good_param_mat = [];
good_data_mat = [];
for tstamp_n = 1:length(t_list)
    curr_tstamp = t_list(tstamp_n);
    curr_param_mati = find(param_mat(:, 1) == curr_tstamp);
    
    %subsets with the same tstamp
    param_mat_tst = param_mat(curr_param_mati, :);
    data_mat_tst = data_mat(curr_param_mati, :);
    
    ar_list = unique(param_mat_tst(:, 2));
    for ar_n = 1:length(ar_list)
        curr_arena = ar_list(ar_n);
        curr_arenai = find(param_mat_tst(:, 2) == curr_arena);
               
        %sub-subsets with the same arena_n
        param_mat_ar = param_mat_tst(curr_arenai, :); 
        data_mat_ar = data_mat_tst(curr_arenai, :); 
        
        cam_list = unique(param_mat_ar(:, 3));
        for cam_n = 1:length(cam_list)
            
        end
        
        
        
    end
    
end

end

function dist_mat = dist_from_center(data_mat, param_mat)
%This function computes distance from center at each time point for each
%fly. dist_mat has n_rows = (n_rows data_mat)/3 (after removing the header
%row) and n_cols = n_cols data_mat (after removing the parameter columns).

n_flies = size(data_mat, 1)./3;
dist_mat = zeros((size(data_mat, 2) - 3),  n_flies) + nan;  %first three columns are parameters, not position measurments        
for fly_n = 1:n_flies
    x_row = ((fly_n - 1).*3) + 1;
    y_row = x_row + 1;
    %making sure to pick right pair of rows
    keyboard
    
    for frame_n = 1:size(data_mat)
        keyboard
        
        
        curr_pair = data_mat(x_row:y_row, frame_n);
        curr_pair = [curr_pair, [0; 0]];
        curr_dist = squareform(pdist(curr_pair'));
        
    end
    
end
end
