clear all
close all

base_paths = [{'C:\Data\Data\Analysed_data\data_sharing\Data_from_Yoshi\PABA_conditioning\MB043C_A1\MB043C'};...
                 {'C:\Data\Data\Analysed_data\data_sharing\Data_from_Yoshi\PABA_conditioning\MB099C_G2\MB099C'}];
             


PE_path = '$20xUAS-CsChrimson-mVenus-attP18$Specificity1-PA-EL-PA+$PA+$rec1\';     %PA v/s EL; PA is CS+
EP_path = '$20xUAS-CsChrimson-mVenus-attP18$Specificity1-PA-EL-PA+$PA+$rec2\';     %PA v/s EL; EL is CS+
PB_path = '$20xUAS-CsChrimson-mVenus-attP18$Specificity3-PA-BA-PA+$PA+$rec1\';     %PA v/s BA; PA is CS+
BP_path = '$20xUAS-CsChrimson-mVenus-attP18$Specificity3-PA-BA-PA+$PA+$rec2\';     %PA v/s BA; BA is CS+
    

frame_time = 1./30; %in s
radius = 91./2;  %in mm

a = colormap('bone'); 
global greymap
greymap = flipud(a);

LED_color = [231, 87, 87]./256;
plus_in_color = LED_color;
plus_out_color = [0.5, 0.9, 0.4];

%User-specified Parameters
t_win = 10;  %trajectory  sampling time window on either side of boundary touch event
sp_thresh = 2;    %in mm/s, the mean speed at which a fly must move in +/-t_win s around a touch event for it to be considered for analysis
rem_win = 5;    %in s, the width of the window to ignore new boundary touch events around one already being counted.
trans_on = 0;   %this boolean switch is to turn translation-alignment of trajectories on or off.
c_dist_thresh = 1;    %threshold fractional distance from center above which quadrant transitions are discarded.
end_t_win = 0.5;  %in s, the time window at the end of the trajectory where mean position is determined to identify returners.

t_win_frs = round(t_win./frame_time);
sp_thresh = sp_thresh.*frame_time;
rem_win_frs = round(rem_win./frame_time);
end_t_win_frs = round(end_t_win./frame_time);

all_plus_in_trajs_E = [];
all_plus_out_trajs_E = [];
all_plus_in_trajs_B = [];
all_plus_out_trajs_B = [];

set_name_vec = [{'PE'}, {'EP'}, {'PB'}, {'BP'}];
for base_path_n = 1:2
    base_path = base_paths{base_path_n, 1};
    if base_path_n == 1
        base_name = 'Alpha1';
    elseif base_path_n == 2
        base_name = 'Gamma2';
    end
    
    if exist([base_path, 'plus_in_traces_aligned_E.mat']) ~= 1
        for set_n = 1:4
            set_name = set_name_vec{1, set_n};
            %identifying current set's CS+ odor
            if strcmp(set_name(1), 'P') == 1        %case where PA was the CS+
                P_plus = 1;
            elseif strcmp(set_name(1), 'P') == 0    %case where PA was not CS+
                P_plus = 0;
            else
            end

            %identifying if current set is fine discrimination or coarse discrimination
            if contains(set_name, 'E') == 0
                E_set = 1;
            elseif contains(set_name, 'E') == 1
                E_set = 0;
            else
            end


            %reading in data
            if exist([base_path, set_name, '_data_struc.mat']) ~= 2
                eval([set_name, '_data_struc = load_data(base_path, ', set_name, '_path, frame_time);']);
                eval(['save([base_path, set_name, ''_data_struc.mat''], ''' set_name '_data_struc'')']);
                eval(['data_struc = ' set_name '_data_struc;']);
            else    
                data_struc = load([base_path, set_name, '_data_struc.mat']);
                eval(['data_struc = data_struc.', set_name, '_data_struc;']);
            end

            %cleaning up data, part1: ensuring consistency, getting rid of spurious
            %points at 0, 0 and adding appropriate offsets to x and y data.
            params_mat = data_struc.params_mat;
            data_mat = data_struc.data_mat;
            
            [params, data_matXY, data_matH, fly_n_replicates, n_disc_rows] = get_good_flies(params_mat, data_mat, base_path_n);

            %reading in or computing distances from center
            if exist([base_path, set_name, '_cdists.mat']) ~= 2
                dist_mat = dist_from_center(data_matXY);
                eval(['save([base_path, ''', set_name, '_cdists.mat''], ''dist_mat'');']);
            else
                dist_mat = load([base_path, set_name, '_cdists.mat']);
                dist_mat = dist_mat.dist_mat;
            end

            %cleaning up data, part2: getting rid of spurious points outside the arena
            [del1, del2] = find(dist_mat > 45.5);
            x_rows = data_matXY(1:2:end, :);
            y_rows = data_matXY(2:2:end, :);
            for pt_n = 1:length(del1)
                x_rows((del1(pt_n)), del2(pt_n)) = nan;
                y_rows((del1(pt_n)), del2(pt_n)) = nan;
            end
            data_matXY = interleave2(x_rows, y_rows, 'row');

            %reading in or sampling trajectories near boundary touches
            if exist([base_path, set_name, '_samp_trajs.mat']) ~= 2
                trajectory_samples = sample_crossing_trajs(data_matXY, t_win_frs, sp_thresh, rem_win_frs, dist_mat, c_dist_thresh, data_struc.od_is(1));
                eval(['save([base_path, ''', set_name, '_samp_trajs.mat''], ''trajectory_samples'');']);
            else
                trajectory_samples = load([base_path, set_name, '_samp_trajs.mat']);
                trajectory_samples = trajectory_samples.trajectory_samples;
            end

            %Pooling and aligning CS+ to CS- transition trajectories and CS- to CS+ transition trajectories 
            [plus_in_trajs, plus_out_trajs] = pool_align_trajs(trajectory_samples, P_plus, (t_win_frs + 1), trans_on);

            if E_set == 1
                %coarse discrimination case, where counter odor was EL
                all_plus_in_trajs_E = [all_plus_in_trajs_E; plus_in_trajs];
                all_plus_out_trajs_E = [all_plus_out_trajs_E; plus_out_trajs];
            elseif E_set == 0
                %fine discrimination case, where counter odor was BA
                all_plus_in_trajs_B = [all_plus_in_trajs_B; plus_in_trajs];
                all_plus_out_trajs_B = [all_plus_out_trajs_B; plus_out_trajs];
            else
            end
        end
        save([base_path, 'plus_in_traces_aligned_E.mat'], 'all_plus_in_trajs_E');
        save([base_path, 'plus_out_traces_aligned_E.mat'], 'all_plus_out_trajs_E');
        save([base_path, 'plus_in_traces_aligned_B.mat'], 'all_plus_in_trajs_B');
        save([base_path, 'plus_out_traces_aligned_B.mat'], 'all_plus_out_trajs_B');
    else
        all_plus_in_trajs_E = load([base_path, 'plus_in_traces_aligned_E.mat']);
        all_plus_in_trajs_E = all_plus_in_trajs_E.all_plus_in_trajs_E;
        all_plus_in_trajs_B = load([base_path, 'plus_in_traces_aligned_B.mat']);
        all_plus_in_trajs_B = all_plus_in_trajs_B.all_plus_in_trajs_B;

        all_plus_out_trajs_E = load([base_path, 'plus_out_traces_aligned_E.mat']);
        all_plus_out_trajs_E = all_plus_out_trajs_E.all_plus_out_trajs_E;
        all_plus_out_trajs_B = load([base_path, 'plus_out_traces_aligned_B.mat']);
        all_plus_out_trajs_B = all_plus_out_trajs_B.all_plus_out_trajs_B;

    end

    %plotting all boundary-touch trajectories
    plot_trajectories(all_plus_out_trajs_E, 'Leaving CS+ quadrant, coarse', 1, 0)
    plot_trajectories(all_plus_in_trajs_E, 'Entering CS+ quadrant, coarse', 1, 0)

    plot_trajectories(all_plus_out_trajs_B, 'Leaving CS+ quadrant, fine', 1, 0)
    plot_trajectories(all_plus_in_trajs_B, 'Entering CS+ quadrant, fine', 1, 0)


    %sub-sampling trajectories that come back and end on the side they started

    %coarse discrimination
    x_rows = all_plus_out_trajs_E(1:3:end, :);
    mean_pos_plus_out = mean(x_rows(:, (end - end_t_win_frs):end), 2, 'omitnan');
    returners_plus_out = find(mean_pos_plus_out < 0);    %identifying flies that end up on the plus side, among plus out trajectories
    returners_plus_out_trajs_E = sub_sample_flies(all_plus_out_trajs_E, returners_plus_out);
    x_rows = all_plus_in_trajs_E(1:3:end, :);
    mean_pos_plus_in = mean(x_rows(:, (end - end_t_win_frs):end), 2, 'omitnan');
    returners_plus_in = find(mean_pos_plus_in > 0);    %identifying flies that end up on the minus side, among plus in trajectories
    returners_plus_in_trajs_E = sub_sample_flies(all_plus_in_trajs_E, returners_plus_in);

    %fine discrimination
    x_rows = all_plus_out_trajs_B(1:3:end, :);
    mean_pos_plus_out = mean(x_rows(:, (end - end_t_win_frs):end), 2, 'omitnan');
    returners_plus_out = find(mean_pos_plus_out < 0);    %identifying flies that end up on the plus side, among plus out trajectories
    returners_plus_out_trajs_B = sub_sample_flies(all_plus_out_trajs_B, returners_plus_out);
    x_rows = all_plus_in_trajs_B(1:3:end, :);
    mean_pos_plus_in = mean(x_rows(:, (end - end_t_win_frs):end), 2, 'omitnan');
    returners_plus_in = find(mean_pos_plus_in > 0);    %identifying flies that end up on the minus side, among plus in trajectories
    returners_plus_in_trajs_B = sub_sample_flies(all_plus_in_trajs_B, returners_plus_in);

    %computing fraction of good returns v/s failed returns for E and B
    %Note that my criterion for returners may not capture all retureners, so
    %these fractions are likely under-representative of the real values

    frac_plus_in_returns_E = (size(returners_plus_in_trajs_E, 1)./3) ./ (size(all_plus_in_trajs_E, 1)./3)
    frac_plus_in_returns_B = (size(returners_plus_in_trajs_B, 1)./3) ./ (size(all_plus_in_trajs_B, 1)./3)

    frac_plus_out_returns_E = (size(returners_plus_out_trajs_E, 1)./3) ./ (size(all_plus_out_trajs_E, 1)./3)
    frac_plus_out_returns_B = (size(returners_plus_out_trajs_B, 1)./3) ./ (size(all_plus_out_trajs_B, 1)./3)


    %plotting returner trajectories
    plot_trajectories(returners_plus_out_trajs_E, ['Leaving CS+ quadrant, coarse ', base_name], 1, 0)
    plot_trajectories(returners_plus_in_trajs_E, ['Entering CS+ quadrant, coarse ', base_name], 1, 0)

    plot_trajectories(returners_plus_out_trajs_B, ['Leaving CS+ quadrant, fine ', base_name], 1, 0)
    plot_trajectories(returners_plus_in_trajs_B, ['Entering CS+ quadrant, fine ', base_name], 1, 0)


    %identifying point of return in returner
    %trajectories. Also point of farthest excursion. Measuring times to these
    %two points.
    %1.farthest excursion
    x_rows = returners_plus_out_trajs_E(1:3:end, t_win_frs:end);
    [max_dists_plus_out_E, fari_plus_out_trajs_E] = max(x_rows, [], 2, 'omitnan');
    max_d_times_plus_out_E = fari_plus_out_trajs_E.*frame_time;     %times in s
    x_rows = returners_plus_in_trajs_E(1:3:end, t_win_frs:end);
    [max_dists_plus_in_E, fari_plus_in_trajs_E] = min(x_rows, [], 2, 'omitnan');
    max_dists_plus_in_E = abs(max_dists_plus_in_E);
    max_d_times_plus_in_E = fari_plus_in_trajs_E.*frame_time;     %times in s

    x_rows = returners_plus_out_trajs_B(1:3:end, t_win_frs:end);
    [max_dists_plus_out_B, fari_plus_out_trajs_B] = max(x_rows, [], 2, 'omitnan');
    max_d_times_plus_out_B = fari_plus_out_trajs_B.*frame_time;     %times in s
    x_rows = returners_plus_in_trajs_B(1:3:end, t_win_frs:end);
    [max_dists_plus_in_B, fari_plus_in_trajs_B] = min(x_rows, [], 2, 'omitnan');
    max_dists_plus_in_B = abs(max_dists_plus_in_B);
    max_d_times_plus_in_B = fari_plus_in_trajs_B.*frame_time;     %times in s


    %2. point of return
    ret_time_in_E = find_ret_frames(returners_plus_in_trajs_E, 'in', frame_time);
    ret_time_in_E = ret_time_in_E.*frame_time;
    ret_time_out_E = find_ret_frames(returners_plus_out_trajs_E, 'out', frame_time);
    ret_time_out_E = ret_time_out_E.*frame_time;

    ret_time_in_B = find_ret_frames(returners_plus_in_trajs_B, 'in', frame_time);
    ret_time_in_B = ret_time_in_B.*frame_time;
    ret_time_out_B = find_ret_frames(returners_plus_out_trajs_B, 'out', frame_time);
    ret_time_out_B = ret_time_out_B.*frame_time;


    %Relating speed to quadrant excursion to identify a transition
    speed_mat = returners_plus_in_trajs_E(3:3:end, :);
    speed_mat = movmean(speed_mat, (0.1./frame_time), 2);


    %Plotting and testing trajectory metrics
    %1. time to return
    data_mat_time = ret_time_in_E;
    data_mat_time = pad_n_concatenate(data_mat_time, ret_time_in_B, 2, nan);
    data_mat_time = pad_n_concatenate(data_mat_time, ret_time_out_E, 2, nan);
    data_mat_time = pad_n_concatenate(data_mat_time, ret_time_out_B, 2, nan);
    fig_h = figure('Name', ['time to return ', base_name]);
    marker_colors = [plus_in_color; plus_in_color; plus_out_color; plus_out_color];
    col_pairs = [];
    line_colors = zeros(4, 3) + 1;    %no lines needed as these are not paired samples
    xlabels = [{'coarse'}, {'fine'}, {'coarse'}, {'fine'}];
    mean_color = [0, 0, 0];
    fig_h = scattered_dot_plot_ttest(data_mat_time, fig_h, 2, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, mean_color, 2, 0.05);
    ylabel('time to return (s)');
    ax_vals = axis;
    ax_vals(1, 4) = 11.5;
    axis(ax_vals);
    fig_wrapup_handle(fig_h, []);

    %2. longest excursion dist.
    data_mat_dist = max_dists_plus_in_E;
    data_mat_dist = pad_n_concatenate(data_mat_dist, max_dists_plus_in_B, 2, nan);
    data_mat_dist = pad_n_concatenate(data_mat_dist, max_dists_plus_out_E, 2, nan);
    data_mat_dist = pad_n_concatenate(data_mat_dist, max_dists_plus_out_B, 2, nan);
    fig_h = figure('Name', ['maximum excursion distance ', base_name]);
    marker_colors = [plus_in_color; plus_in_color; plus_out_color; plus_out_color];
    col_pairs = [];
    line_colors = zeros(4, 3) + 1;    %no lines needed as these are not paired samples
    xlabels = [{'coarse'}, {'fine'}, {'coarse'}, {'fine'}];
    mean_color = [0, 0, 0];
    fig_h = scattered_dot_plot_ttest(data_mat_dist, fig_h, 2, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, mean_color, 2, 0.05);
    ylabel('max excursion dist (mm)');
    ax_vals = axis;
    ax_vals(1, 4) = 47;
    axis(ax_vals);
    fig_wrapup_handle(fig_h, []);



    %3. longest excursion time
    data_mat_time_ex = max_d_times_plus_in_E;
    data_mat_time_ex = pad_n_concatenate(data_mat_time_ex, max_d_times_plus_in_B, 2, nan);
    data_mat_time_ex = pad_n_concatenate(data_mat_time_ex, max_d_times_plus_out_E, 2, nan);
    data_mat_time_ex = pad_n_concatenate(data_mat_time_ex, max_d_times_plus_out_B, 2, nan);
    fig_h = figure('Name',['time to max excursion ', base_name]);
    marker_colors = [plus_in_color; plus_in_color; plus_out_color; plus_out_color];
    col_pairs = [];
    line_colors = zeros(4, 3) + 1;    %no lines needed as these are not paired samples
    xlabels = [{'coarse'}, {'fine'}, {'coarse'}, {'fine'}];
    mean_color = [0, 0, 0];
    fig_h = scattered_dot_plot_ttest(data_mat_time_ex, fig_h, 2, 4, 8, marker_colors, 1, col_pairs, line_colors, xlabels, 1, mean_color, 2, 0.05);
    ylabel('time to max exc. (s)');
    ax_vals = axis;
    ax_vals(1, 4) = 11;
    axis(ax_vals);
    fig_wrapup_handle(fig_h, []);

    
    keyboard

end
%----------------
%worker functions
function [data_struc] = load_data(base_path, spec_path, frame_time)

%reading in and parsing trajectory data from file
[data_mat, params_mat] = read_data_file([base_path, spec_path, 'Results.csv']);

%reformatting data before passing
fnums = data_mat(1, :);
od_oni = find(fnums(1, :) == 0);      %frame_n when LED came on; checked that it's the same for all datasets
od_offi = od_oni + round(60./frame_time);

data_mat(1, :) = [];
data_mat = sortrows(data_mat);
params_mat = data_mat(:, 1:3);        %timestamp, arena_n, cam_n, fly_n and param_n for each row of data_mat
data_mat(:, 1:3) = [];
data_mat = data_mat.*0.1;         %converting pixels into mm (1 px = 0.1mm)

data_struc.data_mat = data_mat;
data_struc.params_mat = params_mat;
data_struc.fnums = fnums;
data_struc.od_is = [od_oni, od_offi];


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
        sub_str = ',';  %curr_line(8);
    else
    end
    
    delimiteri = findstr(curr_line, sub_str);
    curr_row = cell(1, 7);
    %reading in metadata text
    for col_n = 1:7
        curr_delimiters = [(delimiteri(col_n) + 1), (delimiteri(col_n + 1) - 1)];        
        curr_num = str2num(curr_line(curr_delimiters(1):curr_delimiters(2)));
        if isempty(curr_num) == 1
            curr_num = {curr_line(curr_delimiters(1):curr_delimiters(2))};
        else
        end
        curr_row{1, col_n} = curr_num;
    end
    met_data_mat = [met_data_mat; curr_row];    %metadata matrix
    
    %reading in numeric data
    curr_row = zeros(1, (length(delimiteri) - 7)) + nan;
    for col_n = 8:length(delimiteri)
        
        if col_n < length(delimiteri)
            curr_delimiters = [(delimiteri(col_n) + 1), (delimiteri((col_n + 1)) - 1)];
        elseif col_n == length(delimiteri)
            curr_delimiters = [(delimiteri(col_n) + 1), length(curr_line)];
        else
        end
        
        curr_num = str2num(curr_line(curr_delimiters(1):curr_delimiters(2)));
        curr_row(1, (col_n - 7)) = curr_num;
           
    end
   
    %adding fly_n and param_type (x - 0, y - 1, heading - 2) in that order at row beginning to allow easy sorting
    if line_n > 1
        %1. movie_n and fly_n
        col_n = 7;
        curr_delimiters = [(delimiteri(col_n) + 1), (delimiteri(col_n + 1) - 1)];
        curr_str = curr_line(curr_delimiters(1):curr_delimiters(2));
        movie_stri = findstr(curr_str, 'movie');
        fly_stri = findstr(curr_str, 'fly');
        movie_n = str2num(curr_str((movie_stri + 5):(fly_stri - 2)));
        fly_n = str2num(curr_str((fly_stri + 3):end));
        
        %2. param_n
        col_n = 6;
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
                
        params = [movie_n, fly_n, param_n];
        
    elseif line_n == 1
        params = zeros(1, 3) + nan;  %header row
    else
    end
    
    curr_row = [params, curr_row];
    data_mat = [data_mat; curr_row];   
     
end
fclose(fid);
end
end

function [good_param_mat, data_matXY, data_matH, fly_n_replicates, n_disc_rows] = get_good_flies(param_mat, data_mat, path_n)
%This function hierarchically searches param mat for three rows of
%measurements for the same fly, for each fly_n, for each movie_n. It then concatenates such
%good row triplets together to generate a cleaned dataset of observations
%for further analysis.

good_param_mat = [];
good_data_mat = [];
fly_n_replicates = 0;
n_disc_rows = 0;

mov_list = unique(param_mat(:, 1));
for mov_n = 1:length(mov_list)
    curr_mov = mov_list(mov_n);
    curr_movi = find(param_mat(:, 1) == curr_mov);

    %subsets with the same mov_n
    param_mat_mov = param_mat(curr_movi, :);
    data_mat_mov = data_mat(curr_movi, :);

    fly_list = unique(param_mat_mov(:, 2));
    for fly_n = 1:length(fly_list)
        curr_fly = fly_list(fly_n);
        curr_flyi = find(param_mat_mov(:, 2) == curr_fly);

        %sub-sub-sub-subsets with the same fly_n. 
        param_mat_fly = param_mat_mov(curr_flyi, :);
        data_mat_fly = data_mat_mov(curr_flyi, :);

        %This set should contain three rows, one for each measured parameter, else
        %discard this fly.

        if length(unique(param_mat_fly(:, 3))) == 3
            if size(param_mat_fly, 1) > 3       %case where padding rows of zeros exist
                data_summed = sum(data_mat_tr, 2);
                real_rowsi = find(data_summed ~= 0);
                param_mat_fly = param_mat_fly(real_rowsi, :);
                data_mat_fly = data_mat_fly(real_rowsi, :);
            else
            end

            if size(param_mat_fly, 1) > 3
                fly_n_replicates = fly_n_replicates + 1;
                n_disc_rows = n_disc_rows + size(param_mat_tr, 1);
                param_mat_fly = [];
                data_mat_fly = [];
                
            else
            end

            good_param_mat = [good_param_mat; param_mat_fly];
            good_data_mat = [good_data_mat; data_mat_fly];
        else
            %not concatenating and so, discarding current fly
        end


    end
       
    
end
data_mat = good_data_mat;

%getting rid of supurious points that were exactly at 0, 0
x_rows = 1:3:size(data_mat, 1);
y_rows = 2:3:size(data_mat, 1);
h_rows = y_rows + 1;
comb_mat = abs(data_mat(x_rows, :)) + abs(data_mat(y_rows, :));
[del] = find(comb_mat == 0);      %these are fly_ns

x_mat = data_mat(x_rows, :);
x_mat(del) = nan;
y_mat = data_mat(y_rows, :);
y_mat(del) = nan;
data_matXY = interleave2(x_mat, y_mat, 'row');
data_matH = data_mat(h_rows, :).*10;

%correcting offset in co-ordinates
if path_n == 2
    x_rows = 1:2:size(data_matXY, 1);
    y_rows = x_rows + 1;
    data_matXY(x_rows, :) = data_matXY(x_rows, :) - 5;  
    data_matXY(y_rows, :) = data_matXY(y_rows, :) + 5;
else
end
end

function dist_mat = dist_from_center(data_mat)
%This function computes distance from center at each time point for each
%fly. dist_mat has n_rows = (n_rows data_mat)/3 (after removing the header
%row) and n_cols = n_cols data_mat (after removing the parameter columns).

n_flies = size(data_mat, 1)./2;
dist_mat = zeros(n_flies, size(data_mat, 2)) + nan;  %first three columns are parameters, not position measurments        
for fly_n = 1:n_flies
    x_row = ((fly_n - 1).*2) + 1;
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

function [trajectory_samples] = sample_crossing_trajs(data_matXY, t_win_frs, sp_thresh, rem_win_frs, dist_mat, c_dist_thresh, od_oni)
    V_mean_speeds = [];
    H_mean_speeds = [];

    x_rows = data_matXY(1:2:end, od_oni:end);
    y_rows = data_matXY(2:2:end, od_oni:end);
   
    %Identifying quadrant boundary touch points in trajectories (points within 1 mm of quadrant boundaries) 
    %1. Vertical boundaries
    [V_fly_ns, V_frame_ns] = find(abs(x_rows) <0.5);   %vertical boundary transitions
    %classifying boundary touch events and sampling trajectories for +/- t_win seconds
    q12_events = [];    %touch events originating in quadrant1, towards 2 (top left to top right)
    q21_events = [];
    q34_events = [];    %touch events originating in quardant3, towards 4 (bottom right to bottom left)
    q43_events = [];
       
       
    for pt_n = 1:size(V_fly_ns, 1)
        curr_fly_n = V_fly_ns(pt_n);
        curr_fr_n = V_frame_ns(pt_n);
        
        if (V_frame_ns(pt_n) + t_win_frs) > size(x_rows, 2)
            %case where boundary touch happened too close to end of data, discarding touch event.
            continue
        elseif (V_frame_ns(pt_n) - t_win_frs) < 1
            %case where boundary touch happened too close to beginning of data, discarding touch event.
            continue
        elseif dist_mat(curr_fly_n, curr_fr_n) > (45.5).*c_dist_thresh
            %case where boundary touch happened too close to arena edge.
            continue
        else
        end
                
        if isnan(x_rows(V_fly_ns(pt_n), (V_frame_ns(pt_n)))) == 0
           
            x_pts = x_rows(V_fly_ns(pt_n), (V_frame_ns(pt_n) - t_win_frs):(V_frame_ns(pt_n) + t_win_frs) );
            
            x_rows(V_fly_ns(pt_n), V_frame_ns(pt_n):(V_frame_ns(pt_n) + rem_win_frs) ) = nan;          %removing frames rem_win_frs after current event so it's not re-counted
            y_pts = y_rows(V_fly_ns(pt_n), (V_frame_ns(pt_n) - t_win_frs):(V_frame_ns(pt_n) + t_win_frs) );
        else
            continue
        end
        
        %removing trajectories that never cross 0 as ambiguous cross-events
        x_pts_s = x_pts;
        x_pts_s(x_pts_s == 0) = [];
        sign_vec = unique(sign(x_pts_s));
        sign_vec(isnan(sign_vec)) = [];
        if length(sign_vec) < 2
            continue
        else
        end

        %applying a minimum mean speed criterion in a +/-t_win time window around the boundary touch event for further consideration
        %computing speed at each frame (dist from last frame)
        speed_vec = do_speed_check(x_pts, y_pts);
        V_mean_speeds = [V_mean_speeds; mean(speed_vec, 'omitnan')];
        if mean(speed_vec, 'omitnan') < sp_thresh      %in mm/s
            continue
        elseif max(speed_vec, [], 'omitnan') > 2
            continue        %discarding if max speed > 20mm/s - indicative of tracking errors
        else
            %discarding vert-cross trajectories that have post-cross points
            %on both sides of horizontal boundary (ie. multiple cpt crossings)
            unique_sides = unique(sign(y_pts));
            unique_sides(isnan(unique_sides)) = [];
            if length(unique_sides) > 1
                %discarding vert-cross trajectories that have post-cross points
                %on both sides of horizontal boundary (ie. multiple cpt crossings)
                continue
            else
            end
            
        end
        
        %checking if fly started on left or right
        if y_rows(V_fly_ns(pt_n), V_frame_ns(pt_n)) > 0         %case where touch event occured above arena center
            %case where trajectory started on left
            if mean(x_pts((t_win_frs - 10):t_win_frs)) < 0      %computing mean position in 1s before touch event     
                q12_events = [q12_events; x_pts; y_pts; speed_vec];
            elseif mean(x_pts((t_win_frs - 10):t_win_frs)) > 0    %case where trajectory started on the right
                q21_events = [q21_events; x_pts; y_pts; speed_vec];
            else
            end

        elseif y_rows(V_fly_ns(pt_n), V_frame_ns(pt_n)) < 0     %case where touch event occured below arena center
            %case where trajectory started on left
            if mean(x_pts((t_win_frs - 10):t_win_frs)) < 0      %computing mean position in 1s before touch event     
                q43_events = [q43_events; x_pts; y_pts; speed_vec];
            elseif mean(x_pts((t_win_frs - 10):t_win_frs)) > 0    %case where trajectory started on the right
                q34_events = [q34_events; x_pts; y_pts; speed_vec];
            else
            end
        end

    end
    
    x_rows = data_matXY(1:2:end, :);
    y_rows = data_matXY(2:2:end, :);
    y_rows_mk = y_rows;
    [H_fly_ns, H_frame_ns] = find(abs(y_rows) <0.5);   %horixontal boundary transitions
    %classifying boundary touch events and sampling trajectories for +/- t_win seconds
    curr_fly_n = V_fly_ns(pt_n);
    curr_fr_n = V_frame_ns(pt_n);
    q14_events = [];    %touch events originating in quadrant1, towards 2 (top left to top right)
    q41_events = [];
    q23_events = [];    %touch events originating in quardant3, towards 4 (bottom right to bottom left)
    q32_events = [];
    for pt_n = 1:size(H_fly_ns, 1)
        if (H_frame_ns(pt_n) + t_win_frs) > size(x_rows, 2)
            %case where boundary touch happened too close to end of data, discarding touch event.
            continue
        elseif (H_frame_ns(pt_n) - t_win_frs) < 1
            %case where boundary touch happened too close to beginning of data, discarding touch event.
            continue
            elseif dist_mat(curr_fly_n, curr_fr_n) > (45.5).*c_dist_thresh
            %case where boundary touch happened too close to arena edge.
            continue
        else
        end
        
        if isnan(y_rows_mk(H_fly_ns(pt_n), (H_frame_ns(pt_n)))) == 0
            x_pts = x_rows(H_fly_ns(pt_n), (H_frame_ns(pt_n) - t_win_frs):(H_frame_ns(pt_n) + t_win_frs) );
            y_rows_mk(H_fly_ns(pt_n), H_frame_ns(pt_n):(H_frame_ns(pt_n) + rem_win_frs) ) = nan;          %removing current event so it's not re-counted
            y_pts = y_rows(H_fly_ns(pt_n), (H_frame_ns(pt_n) - t_win_frs):(H_frame_ns(pt_n) + t_win_frs) );
        else
            continue
        end
        %removing trajectories that never cross 0 as ambiguous cross-events
        y_pts_s = x_pts;
        y_pts_s(y_pts_s == 0) = [];
        sign_vec = unique(sign(y_pts_s));
        sign_vec(isnan(sign_vec)) = [];
        if length(sign_vec) < 2
            continue
        else
        end

        
        
        %applying a minimum mean speed criterion in a +/-t_win time window around the boundary touch event for further consideration
        %computing speed at each frame (dist from last frame)
        speed_vec = do_speed_check(x_pts, y_pts);
        H_mean_speeds = [H_mean_speeds; mean(speed_vec, 'omitnan')];
        if mean(speed_vec, 'omitnan') < sp_thresh      %in mm/s
            continue
        elseif max(speed_vec, [], 'omitnan') > 2
            continue        %discarding if max speed > 20mm/s - indicative of tracking errors
        else
            unique_sides = unique(sign(x_pts));
            unique_sides(isnan(unique_sides)) = [];
            if length(unique_sides) > 1
                %discarding vert-cross trajectories that have post-cross points
                %on both sides of horizontal boundary (ie. multiple cpt crossings)
                continue
            else
            end
        end
        
        %checking if fly started on top or bottom
        if x_rows(H_fly_ns(pt_n), H_frame_ns(pt_n)) < 0     %case where boundary touch happened right of center
            %case where trajectory started on top
            if mean(y_pts((t_win_frs - 10):t_win_frs)) > 0      %computing mean position in 1s before touch event     
                q14_events = [q14_events; x_pts; y_pts; speed_vec];
            elseif mean(y_pts((t_win_frs - 10):t_win_frs)) < 0    %case where trajectory started on the right
                q41_events = [q41_events; x_pts; y_pts; speed_vec];
            else
            end

        elseif x_rows(H_fly_ns(pt_n), H_frame_ns(pt_n)) > 0 %case where boundary touch happened left of center
            %case where trajectory started on top
            if mean(y_pts((t_win_frs - 10):t_win_frs)) > 0      %computing mean position in 1s before touch event     
                q23_events = [q23_events; x_pts; y_pts; speed_vec];
            elseif mean(y_pts((t_win_frs - 10):t_win_frs)) < 0    %case where trajectory started on the bottom
                q32_events = [q32_events; x_pts; y_pts; speed_vec];
            else
            end
        end
    end
    trajectory_samples.q12 = q12_events;
    trajectory_samples.q21 = q21_events;
    trajectory_samples.q34 = q34_events;
    trajectory_samples.q43 = q43_events;
    
    trajectory_samples.q14 = q14_events;
    trajectory_samples.q41 = q41_events;
    trajectory_samples.q23 = q23_events;
    trajectory_samples.q32 = q32_events;

function [speed_vec] = do_speed_check(x_pts, y_pts)
    speed_mat = pdist([x_pts', y_pts']);
    speed_mat = squareform(speed_mat);
    speed_vec = zeros(1, size(x_pts, 2)) + nan;
    for pt_ni = 1:(size(x_pts, 2) - 1)
        speed_vec(1, (pt_ni + 1)) = speed_mat(pt_ni, (pt_ni + 1));
    end
    
end
   
end

function speed_mat = calc_speed(data_mat)
    %smoothing position traces to avoid jumpy speed measurments
    data_mat = movmean(data_mat', 5)';
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

function [plus_in_trajs, plus_out_trajs] = pool_align_trajs(trajectory_samples, P_plus, ref_fr, trans_on)
    P_out_trajs = [];
    P_in_trajs = [];
    %applying transformations to align 
    %assuming PA is CS+ and correcting appropriately at the end
    q12_trajs = trajectory_samples.q12;
    trajs_adj = apply_rot_trans(q12_trajs, 0, 0, ref_fr, trans_on);
    P_out_trajs = [P_out_trajs; trajs_adj];
    
    q14_trajs = trajectory_samples.q14;
    trajs_adj = apply_rot_trans(q14_trajs, 90, 1, ref_fr, trans_on);
    P_out_trajs = [P_out_trajs; trajs_adj];
    
    q32_trajs = trajectory_samples.q32;
    trajs_adj = apply_rot_trans(q32_trajs, -90, 1, ref_fr, trans_on);
    P_out_trajs = [P_out_trajs; trajs_adj];
    
    q34_trajs = trajectory_samples.q34;
    trajs_adj = apply_rot_trans(q34_trajs, -180, 0, ref_fr, trans_on);
    P_out_trajs = [P_out_trajs; trajs_adj];

    
    q21_trajs = trajectory_samples.q21;
    trajs_adj = apply_rot_trans(q21_trajs, 0, 0, ref_fr, trans_on);
    P_in_trajs = [P_in_trajs; trajs_adj];
     
    q23_trajs = trajectory_samples.q23;
    trajs_adj = apply_rot_trans(q23_trajs, -90, 1, ref_fr, trans_on);
    P_in_trajs = [P_in_trajs; trajs_adj];
    
    q41_trajs = trajectory_samples.q41;
    trajs_adj = apply_rot_trans(q41_trajs, 90, 1, ref_fr, trans_on);
    P_in_trajs = [P_in_trajs; trajs_adj];
    
    q43_trajs = trajectory_samples.q43;
    trajs_adj = apply_rot_trans(q43_trajs, -180, 0, ref_fr, trans_on);
    P_in_trajs = [P_in_trajs; trajs_adj];
    
    %flipping y co-ords if PA was CS- to ensure CS + is the quadrant above
    if P_plus == 0
        mult_mat = ones(size(P_in_trajs, 1), size(P_in_trajs, 2));
        mult_mat(1:3:end, :) = -1;
        P_in_trajs = P_in_trajs.*mult_mat;
        
        mult_mat = ones(size(P_out_trajs, 1), size(P_out_trajs, 2));
        mult_mat(1:3:end, :) = -1;
        P_out_trajs = P_out_trajs.*mult_mat;
        
        plus_in_trajs = P_out_trajs;
        plus_out_trajs = P_in_trajs;
    else
        plus_in_trajs = P_in_trajs;
        plus_out_trajs = P_out_trajs;
    end
    
%     %testing lines
%     figure(1)
%     x_rows = P_in_trajs(1:3:end, :);
%     y_rows = P_in_trajs(2:3:end, :);
%     plot(x_rows(:, 1:ref_fr)', y_rows(:, 1:ref_fr)', 'Color', [0.65, 0.65, 0.65]);
%     hold on
%     plot(x_rows(:, ref_fr:end)', y_rows(:, ref_fr:end)', 'Color', 'k');
%     hold off
%     
%     figure(2)
%     x_rows = P_out_trajs(1:3:end, :);
%     y_rows = P_out_trajs(2:3:end, :);
%     plot(x_rows(:, 1:ref_fr)', y_rows(:, 1:ref_fr)', 'Color', [0.65, 0.65, 0.65]);
%     hold on
%     plot(x_rows(:, ref_fr:end)', y_rows(:, ref_fr:end)', 'Color', 'k');
%     hold off
%     keyboard
% 
    
end

function [trajs_out] = apply_rot_trans(trajs_in, rotation, flip, ref_fr, trans_on)
    %converts x, y co-ordinate rows, to polar coordinates,
    %rotates by angle rotation, converts back to cartesian coords and flips 
    %vertically if flip == 1. Then subtracts x, y co-ords at specified
    %frame_n to translate and align.
    trajs_out = [];
    x_rows = trajs_in(1:3:end, :);
    y_rows = trajs_in(2:3:end, :);
    sp_rows = trajs_in(3:3:end, :);
   
    [theta, rho] = cart2pol(x_rows, y_rows);
    
    rot_ang = deg2rad(rotation);
    theta = theta - rot_ang;
    [x_rows, y_rows] = pol2cart(theta, rho);
    
    if trans_on == 1
        %subtracting x, y coords at reference frame
        ref_fr_x = x_rows(:, ref_fr);
        x_rows_adj = x_rows - repmat(ref_fr_x, 1, size(x_rows, 2)) ;
        ref_fr_y = y_rows(:, ref_fr);
        y_rows_adj = y_rows - repmat(ref_fr_y, 1, size(y_rows, 2));
    else
        x_rows_adj = x_rows;
        y_rows_adj = y_rows;
    end
    
    if flip == 1
        x_rows_adj = x_rows_adj.*(-1);
    else
    end
    
    trajs_out = interleave2(x_rows_adj, y_rows_adj, sp_rows, 'row');
    
end

function [] = plot_trajectories(traj_mat, name, col_type, zoomin)

    if col_type == 1
        col1 = [0.65, 0.65, 0.65];
        col2 = [0, 0, 0];
    else
    end

    fig_h = figure('Name', name);
    mid_fr = ceil(size(traj_mat, 2)./2);
    x_rows = traj_mat(1:3:end, :);
    y_rows = traj_mat(2:3:end, :);
    plot(x_rows(:, 1:(mid_fr -1))', y_rows(:, 1:(mid_fr - 1))', 'Color', col1)
    hold on
    plot(x_rows(:, mid_fr:end)', y_rows(:, mid_fr:end)', 'Color', col2)
    plot([1, 1], [-50, 50], '--', 'Color', [231, 87, 87]./256)
    plot([-1, -1], [-50, 50], '--', 'Color', [231, 87, 87]./256)
    
    if zoomin == 0
        axis([-50, 50, 0, 50])
    elseif zoomin == 1
        axis([-10, 10, -10, 10])
    else
    end
    
    ylabel('distance (mm)');
    xlabel('distance (mm)');
    pbaspect([1, 0.5, 1]);
    hold off
    fig_wrapup_handle(fig_h, []);
    text(-40, 40, 'CS+');
    text(30, 40, 'CS-');

end

function [sub_traj_mat] = sub_sample_flies(traj_mat, fly_n_vec)
    
    n_flies = length(fly_n_vec);
    sub_traj_mat = zeros((n_flies.*3), size(traj_mat, 2)) + nan;
    for fly_n = 1:n_flies
        fly_ni = fly_n_vec(fly_n);
        curr_row = (fly_ni - 1).*3 + 1;
        curr_rows_from = [curr_row; (curr_row + 1); (curr_row + 2)]; %row numbers in source mat for x, y and speed for current fly
        
        curr_row = (fly_n - 1).*3 + 1;
        curr_rows_to = [curr_row; (curr_row + 1); (curr_row + 2)]; %row numbers in sink mat for x, y and speed for current fly
        sub_traj_mat(curr_rows_to, :) = traj_mat(curr_rows_from, :);
        
    end

end

function [ret_frame_n] = find_ret_frames(traj_mat, CS_pl_dir, frame_time)
mid_fr = ceil(size(traj_mat, 2)./2);
x_rows = traj_mat(1:3:end, mid_fr:end);     %mid_fr is the first cross-over frame. Looking for return after that.
x_rows = movmean(x_rows, round(0.25./frame_time), 2);
if strcmp(CS_pl_dir, 'in') == 0
    x_rows = x_rows.*(-1);
else
end
x_rows_signi = sign(x_rows);            %binarising by side of zero
x_rows_difi = diff(x_rows_signi, 1, 2); %looking for changes in side of zero
x_rows_difi = x_rows_difi(:, 3:end);    %discarding points too close to beginning
[rfly_n, ret_frame_n] = find(x_rows_difi == 2);
if length(rfly_n) > size(x_rows_difi, 1)
    [del, first_ret, del2] = unique(rfly_n);
    ret_frame_n = ret_frame_n(first_ret);
else
end
end
