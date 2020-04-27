clear all
close all

base_path = 'C:\Data\Data\Analysed_data\data_sharing\Data_from_Yoshi\upwind_screen_trajectories\';
fname_empty = 'Results_10sLED-air_Empty-split';
fname_MB077B = 'Results_10sLED-air_MB077B';

frame_time = 1./30; %in s
dia = 910;  %in pixels

a = colormap('bone');
global greymap
greymap = flipud(a);

%reading in parsed data
[ctrl_data_mat, ctrl_params, ctrl_fnums, MBON_data_mat, MBON_params, MBON_fnums, LED_oni, LED_offi] = load_data(base_path, fname_empty, fname_MB077B, frame_time);

%cleaning up data
[ctrl_params, ctrl_data_mat, ctrl_fly_n_replicates, crtl_n_disc_rows] = get_good_flies(ctrl_params, ctrl_data_mat);
[MBON_params, MBON_data_mat, MBON_fly_n_replicates, MBON_n_disc_rows] = get_good_flies(MBON_params, MBON_data_mat);
fig_n = 0;
close all

fig_n = fig_n + 1;
figure('Name', 'Ex trajectories, MB077B')
t_win = [(LED_offi - round(4./frame_time)), (LED_offi + round(4./frame_time))];
adj_data_mat = plot_trajectories(fig_n, MBON_data_mat, t_win, 0.05, 2);
keyboard
fig_n = fig_n + 1;
figure('Name', 'Ex trajectories, control')
t_win = [LED_offi, (LED_offi + round(8./frame_time))];
adj_data_mat = plot_trajectories(fig_n, ctrl_data_mat, t_win, 0.05, 2);
keyboard

% %getting distances from arena center
% ctrl_c_dists = dist_from_center(ctrl_data_mat);
% ctrl_c_dists = ctrl_c_dists./(dia./2);  %normalising to radius
% stim_frs = [LED_oni, LED_offi];
% 
% %ctrl raw distance to center
% fig_n = fig_n + 1;
% figure('Name', 'control trajectories, d_c_e_n_t_r_e')
% imagesc(ctrl_c_dists)
% colormap(greymap);
% set_xlabels_time(fig_n, frame_time, 10)
% ylabel('fly n')
% fig_wrapup(fig_n, [])
% add_stim_bar(fig_n, stim_frs, [0.8, 0.4, 0.5]);
% 
% %ctrl mean raw distance to center
% fig_n = fig_n + 1;
% figure('Name', 'control mean+/-SE control trajectory')
% ctrl_mean_traj = mean(ctrl_c_dists, 1);
% ctrl_se_traj = std(ctrl_c_dists)./sqrt(size(ctrl_mean_traj, 1));
% shadedErrorBar([], ctrl_mean_traj, ctrl_se_traj, {'Color', [0.65, 0.65, 0.65]})
% set_xlabels_time(fig_n, frame_time, 10)
% ylabel('norm., baseline-subtracted d_c_e_n_t_r_e')
% fig_wrapup(fig_n, [])
% add_stim_bar(fig_n, stim_frs, [0.8, 0.4, 0.5]);
% 
% 
% MBON_c_dists = dist_from_center(MBON_data_mat);
% MBON_c_dists = MBON_c_dists./(dia./2);  %normalising to radius
% 
% MBON_baselines = mean(MBON_c_dists(:, (LED_oni - round(10./frame_time)):(LED_oni - 1)), 2);
% MBON_c_dists_sub  = MBON_c_dists - repmat(MBON_baselines, 1, size(MBON_c_dists, 2));
% 
% %MBON raw distance to center
% fig_n = fig_n + 1;
% figure('Name', 'MB077B trajectories, d_c_e_n_t_r_e')
% imagesc(MBON_c_dists)
% colormap(greymap);
% set_xlabels_time(fig_n, frame_time, 10)
% ylabel('fly n')
% fig_wrapup(fig_n, [])
% add_stim_bar(fig_n, stim_frs, [0.8, 0.4, 0.5]);
% 
% %MBON mean raw distance to center
% fig_n = fig_n + 1;
% figure('Name','MB077B mean+/-SE control trajectory')
% MBON_mean_traj = mean(MBON_c_dists, 1);
% MBON_se_traj = std(MBON_c_dists)./sqrt(size(MBON_mean_traj, 1));
% shadedErrorBar([], MBON_mean_traj, MBON_se_traj, {'Color', [0.65, 0.65, 0.65]})
% set_xlabels_time(fig_n, frame_time, 10)
% ylabel('norm., baseline-subtracted d_c_e_n_t_r_e')
% fig_wrapup(fig_n, []);
% add_stim_bar(fig_n, stim_frs, [0.8, 0.4, 0.5]);

% %ctrl heading
% fig_n = fig_n + 1;
% figure('Name','Ctrl heading')
% ctrl_heading_mat = ctrl_data_mat(3:3:end, :);
% imagesc(ctrl_heading_mat);
% colormap(greymap);
% set_xlabels_time(fig_n, frame_time, 10)
% ylabel('fly n')
% fig_wrapup(fig_n, [])
% add_stim_bar(fig_n, stim_frs, [0.8, 0.4, 0.5]);
% 
% %control mean heading
% fig_n = fig_n + 1;
% figure('Name','control mean+/-SE heading')
% ctrl_mean_head = mean(ctrl_data_mat(3:3:end, :), 1);
% ctrl_se_head = std(ctrl_data_mat(3:3:end, :))./sqrt(size(ctrl_mean_head, 1));
% shadedErrorBar([], ctrl_mean_head, ctrl_se_head, {'Color', [0.65, 0.65, 0.65]})
% set_xlabels_time(fig_n, frame_time, 10)
% ylabel('heading')
% fig_wrapup(fig_n, []);
% add_stim_bar(fig_n, stim_frs, [0.8, 0.4, 0.5]);
% 
% %MBON heading
% fig_n = fig_n + 1;
% figure('Name','MB077B heading')
% MBON_heading_mat = MBON_data_mat(3:3:end, :);
% imagesc(MBON_heading_mat);
% colormap(greymap);
% set_xlabels_time(fig_n, frame_time, 10)
% ylabel('fly n')
% fig_wrapup(fig_n, [])
% add_stim_bar(fig_n, stim_frs, [0.8, 0.4, 0.5]);

% %MBON mean heading
% fig_n = fig_n + 1;
% figure('Name','MB077B mean+/-SE Heading angle')
% MBON_mean_angle = mean(MBON_heading_mat, 1);
% MBON_se_angle = std(MBON_heading_mat)./sqrt(size(MBON_mean_angle, 1));
% shadedErrorBar([], MBON_mean_angle, MBON_se_angle, {'Color', [0.65, 0.65, 0.65]})
% set_xlabels_time(fig_n, frame_time, 10)
% ylabel('norm., baseline-subtracted angle')
% fig_wrapup(fig_n, []);
% add_stim_bar(fig_n, stim_frs, [0.8, 0.4, 0.5]);

%Plotting fly speeds
ctrl_speed_mat = calc_speed(ctrl_data_mat);
MBON_speed_mat = calc_speed(MBON_data_mat);

%ctrl all traces
fig_n = fig_n + 1;
figure('Name','ctrl speed')
imagesc(ctrl_speed_mat);
colormap(greymap);
set_xlabels_time(fig_n, frame_time, 10)
ylabel('fly n')
fig_wrapup(fig_n, [])
add_stim_bar(fig_n, stim_frs, [0.8, 0.4, 0.5]);

%ctrl mean trace
fig_n = fig_n + 1;
figure('Name','ctrl mean+/-SE speed')
ctrl_mean_speed = mean(ctrl_speed_mat, 1);
ctrl_se_speed = std(ctrl_speed_mat)./sqrt(size(ctrl_mean_speed, 1));
shadedErrorBar([], ctrl_mean_speed, ctrl_se_speed, {'Color', [0.65, 0.65, 0.65]})
set_xlabels_time(fig_n, frame_time, 10)
ylabel('norm., baseline-subtracted angle')
fig_wrapup(fig_n, []);
add_stim_bar(fig_n, stim_frs, [0.8, 0.4, 0.5]);

%MBON all traces
fig_n = fig_n + 1;
figure('Name','MB077B speed')
imagesc(MBON_speed_mat);
colormap(greymap);
set_xlabels_time(fig_n, frame_time, 10)
ylabel('fly n')
fig_wrapup(fig_n, [])
add_stim_bar(fig_n, stim_frs, [0.8, 0.4, 0.5]);

%MBON mean trace
fig_n = fig_n + 1;
figure('Name','MB077B mean+/-SE speed')
MBON_mean_speed = mean(MBON_speed_mat, 1);
MBON_se_speed = std(MBON_speed_mat)./sqrt(size(MBON_mean_speed, 1));
shadedErrorBar([], MBON_mean_speed, MBON_se_speed, {'Color', [0.65, 0.65, 0.65]})
set_xlabels_time(fig_n, frame_time, 10)
ylabel('norm., baseline-subtracted angle')
fig_wrapup(fig_n, []);
add_stim_bar(fig_n, stim_frs, [0.8, 0.4, 0.5]);


%----------------
%worker functions
function [ctrl_data_mat, ctrl_params, ctrl_fnums, MBON_data_mat, MBON_params, MBON_fnums, LED_oni, LED_offi] = load_data(base_path, fname_empty, fname_MB077B, frame_time)

if exist([base_path, fname_empty, '_data.mat']) ~= 2
    
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
ctrl_params = ctrl_data_mat(:, 1:6);        %timestamp, arena_n, cam_n, fly_n and param_n for each row of data_mat
MBON_params = MBON_data_mat(:, 1:6);        %timestamp, arena_n, cam_n, fly_n and param_n for each row of data_mat
ctrl_data_mat(:, 1:6) = [];
MBON_data_mat(:, 1:6) = [];

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
        
        %6. trial_n
        col_n = 6;
        curr_delimiters = [(delimiteri(col_n) + 7), (delimiteri(col_n + 1) - 1)];
        trial_n = str2num(curr_line(curr_delimiters(1):curr_delimiters(2)));
        
        params = [tstamp, arena_n, cam_n, fly_n, trial_n, param_n];
    elseif line_n == 1
        params = zeros(1, 6) + nan;  %header row
    else
    end
    curr_row = [params, curr_row];
    data_mat = [data_mat; curr_row];
       
end
fclose(fid);
end
end

function [good_param_mat, good_data_mat, fly_n_replicates, n_disc_rows] = get_good_flies(param_mat, data_mat)
%This function hierarchically searches param mat for three rows of
%measurements for the same fly, for each fly. It then concatenates such
%good row triplets together to generate a cleaned dataset of observations
%for further analysis.

t_list = unique(param_mat(:, 1));   %list of unique time stamps
good_param_mat = [];
good_data_mat = [];
fly_n_replicates = 0;
n_disc_rows = 0;
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
            curr_cam = cam_list(cam_n);
            curr_cami = find(param_mat_ar(:, 3) == curr_cam);
            
            %sub-sub-subsets with the same cam_n
            param_mat_cam = param_mat_ar(curr_cami, :);
            data_mat_cam = data_mat_ar(curr_cami, :);
            
            fly_list = unique(param_mat_cam(:, 4));
            for fly_n = 1:length(fly_list)
                curr_fly = fly_list(fly_n);
                curr_flyi = find(param_mat_cam(:, 4) == curr_fly);
                
                %sub-sub-sub-subsets with the same fly_n. 
                param_mat_fly = param_mat_cam(curr_flyi, :);
                data_mat_fly = data_mat_cam(curr_flyi, :);
                
                tr_list = unique(param_mat_fly(:, 5));
                for trial_n = 1:length(tr_list)
                    curr_tr = tr_list(trial_n);
                    curr_tri = find(param_mat_fly(:, 5) == curr_tr);
                    
                    %sub-sub-sub-subsetsa with the same trial_n.
                    param_mat_tr = param_mat_fly(curr_tri, :);
                    data_mat_tr = data_mat_fly(curr_tri, :);
                    
                    %This set should contain three rows, one for each measured parameter, else
                    %discard this fly.

                    if length(unique(param_mat_tr(:, 6))) == 3
                        if size(param_mat_tr, 1) > 3       %case where padding rows of zeros exist
                            data_summed = sum(data_mat_tr, 2);
                            real_rowsi = find(data_summed ~= 0);
                            param_mat_tr = param_mat_tr(real_rowsi, :);
                            data_mat_tr = data_mat_tr(real_rowsi, :);
                        else
                        end

                        if size(param_mat_tr, 1) > 3
                            fly_n_replicates = fly_n_replicates + 1;
                            n_disc_rows = n_disc_rows + size(param_mat_tr, 1);
                            param_mat_tr = [];
                            data_mat_tr = [];
                        else
                        end

                        good_param_mat = [good_param_mat; param_mat_tr];
                        good_data_mat = [good_data_mat; data_mat_tr];
                    else
                        %not concatenating and so, discarding current fly
                    end
                    
                end
            end
        end
        
        
        
    end
    
end

end

function dist_mat = dist_from_center(data_mat)
%This function computes distance from center at each time point for each
%fly. dist_mat has n_rows = (n_rows data_mat)/3 (after removing the header
%row) and n_cols = n_cols data_mat (after removing the parameter columns).

n_flies = size(data_mat, 1)./3;
dist_mat = zeros(n_flies, size(data_mat, 2)) + nan;  %first three columns are parameters, not position measurments        
for fly_n = 1:n_flies
    x_row = ((fly_n - 1).*3) + 1;
    y_row = x_row + 1;
    
    for frame_n = 1:size(data_mat, 2)
        curr_pair = data_mat(x_row:y_row, frame_n);
        curr_pair = [curr_pair, [0; 0]];
        curr_dist = squareform(pdist(curr_pair'));
        curr_dist = curr_dist(1, 2);
        
        dist_mat(fly_n, frame_n) = curr_dist;
    end
    
end
end

function speed_mat = calc_speed(data_mat)
    n_flies = size(data_mat, 1)./3;
    speed_mat = [];
    n_frames = size(data_mat, 2);
    for fly_n = 1:n_flies
        x_row = ((fly_n - 1).*3) + 1;
        y_row = x_row + 1;
        for frame_n = 1:(n_frames - 1)
            pair1 = data_mat(x_row:y_row, frame_n);
            pair2 = data_mat(x_row:y_row, (frame_n + 1) );
            pairs = [pair1, pair2];
            curr_dist = squareform(pdist(pairs'));
            curr_dist = curr_dist(1, 2);
            speed_mat(fly_n, frame_n) = curr_dist;            
        end        
    end
end

function [adj_data_mat] = plot_trajectories(fig_n, data_mat, t_win, sparseness, col_type)
    adj_data_mat = zeros((size(data_mat, 1).*(2./3)), (t_win(2) - t_win(1) + 1)) + nan;
    %converting xy coords to rtheta coords
    for fly_n = 1:(size(data_mat, 1)./3)
        x_row = ((fly_n - 1).*3 + 1);
        y_row = x_row + 1;
        x_rowi = (fly_n - 1).*2 + 1;
        y_rowi = x_rowi + 1;
        curr_xy = data_mat(x_row:y_row, t_win(1):t_win(2));
        [theta, rho] = cart2pol(curr_xy(1, :), curr_xy(2, :));
        %subtracting initial heading and position to line up all
        %trajetories with each other
        theta = theta - theta(1);
        rho = rho - rho(1);
        [adj_x, adj_y] = pol2cart(theta, rho);
        adj_data_mat(x_rowi:y_rowi, :) = [adj_x; adj_y];
    end
    
            
    %plotting
    %making color vector
    r_vals = linspace(0.9, 0.1, size(adj_data_mat, 2))';
    g_vals = zeros(size(r_vals, 1), 1) + 0.2;
    b_vals = linspace(0.1, 0.9, size(adj_data_mat, 2))';
    col_vec_graded = [r_vals, g_vals, b_vals];
    half_l = round(size(adj_data_mat, 2)./2);
    col_vec_dig = [repmat([0.9, 0.2, 0.1], half_l, 1); repmat([0.1, 0.2, 0.9], half_l, 1)];
    if col_type == 1
        col_vec = col_vec_graded;
    elseif col_type == 2
        col_vec = col_vec_dig;
    else
    end
    figure(fig_n);
    for subplot_n = 1:3
        subplot_tight(1, 3, subplot_n)
        n_flies_r = round(size(data_mat, 1).*sparseness);
        rng('shuffle');
        r_vec = randperm((size(data_mat, 1)./3), n_flies_r);
        
        for r_fly_n = 1:length(r_vec)
            curr_fly = r_vec(r_fly_n);
            x_row = (curr_fly - 1).*2 + 1;
            y_row = x_row + 1;
            for frame_n = 1:size(adj_data_mat, 2)
                curr_color = col_vec(frame_n, :);
                plot(adj_data_mat(x_row, frame_n), adj_data_mat(y_row, frame_n), '.', 'Color', curr_color);
                hold on
            end
            
%             if rem(r_fly_n, 5) == 1
%                 keyboard
%             else
%             end
        end
        hold off
        
    end
        
end
