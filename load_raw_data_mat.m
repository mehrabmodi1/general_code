function [raw_data_mat direc set_type rand_times_list trial_time no_frames no_cells no_trials frame_time...
    CS_onset_time CS_onset_frame CS_duration CS_US_delay US_onset_frame US_duration trial_type_vec learned blink_list chance_blinks ROI_mat] = load_raw_data_mat(list_direc, an_direc, data_direc, rho_control)
%syntax: [] = load_raw_data_mat(direc) 
%load_raw_data_mat loads the data extracted from images and stored in
%raw_data_mat.mat into raw_data_mat. It also sets up all other
%initialisation variables required for further analysis. The inputs it
%requires are list_direc, which is the directory containing the original
%image files from which raw_data_mat.mat was originally extracted, and
%an_direc, which is the generic, analysis directory path within which all
%analysis data is saved (to the level of the date, not inclusive).

if nargin == 1
    an_direc = 'C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\';
    data_direc = 'C:\Data\data\BlinkView\';
    rho_control = 1;
elseif nargin == 2
    data_direc = 'C:\Data\data\BlinkView\';
    rho_control = 1;
elseif nargin == 3
    rho_control = 1;
else
end

mag_an_direc = 'C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\Magnet\blink_magnet_trials_classified\';

datei = findstr(list_direc, '201');
datei = datei(1,1);
    
%checking for the existence of an assembled direc
if isdir([list_direc(1, 1:(end-4)) '_assembled']) == 0
    direc = [an_direc list_direc(datei:end) '\'];
elseif isdir([list_direc(1, 1:(end-4)) '_assembled']) == 1
    direc = [an_direc list_direc(datei:end-4) '_assembled\'];
end
clear datei

%checking if direc is a pseudo-rand dataset
randi = findstr(direc, 'rand');
if isempty(randi) == 1
    set_type = 0;               %0 for trace dataset
elseif isempty(randi) == 0
    set_type = 1;               %1 for rand dataset
end

clear randi

%loading list of trial types and info
%trial_type_vec = load([direc 'trial_type_list.txt']); 
%Temporarily specified below, as trial_type vectors haven't been saved
%correctly by in_vivo_dff_trace_maker
info = load([direc 'info.mat']);
info = info.info;

%parsing direc to find dataset name
datei = findstr(direc, '201');
datei = datei(1,1);
set_name = direc(1, datei:end);
uscorei = findstr(set_name, '_');
control_name = [set_name(1, 1:uscorei) 'control_no-puff-ROI\'];
mag_set_name = [direc(1, datei:(datei + 7) ), direc(1, (datei + 9):end )];

if rho_control == 0
    set_name = control_name;
    direc = ['C:\Stuff\NCBS\Lab_13\Data\analysis\Analysed Data\BlinkView\' set_name];        
elseif rho_control == 1
end

clear uscorei


%reading in file with extracted image data
raw_data_mat = load([direc 'raw_int_data.mat']);
raw_data_mat = raw_data_mat.raw_data_mat;

ROI_mat = load([direc 'ROIs.txt']);         %reading in original ROIs used to extract data

%reading in list of random times if dataset is pseudo-rand
if set_type == 1
    rand_times_list = load([data_direc set_name 'rand CS timings.xls']);
else
    rand_times_list = [];
end

clear datei

%adjusting for old format of info.txt
if info(1, 1) == 15 | info(1, 1) == 25 | info(1, 1) == 6
    info = info(2:end, 1);
else
end

%finding out length of each trial in s
temp = load([data_direc set_name 'Trial - 0-1.xls']);
trial_time = size(temp, 1)./10000;              %since acquisition rate of blink trace is 10KHz
clear temp

no_frames = info(1, 1);
frame_time = (trial_time./no_frames)*1000;          %in ms

no_trials = size(raw_data_mat, 3);
trial_type_vec = zeros(no_trials, 1);
no_cells = size(raw_data_mat, 2);
CS_onset_time = info(4, 1);                         %in ms
CS_onset_frame = round(CS_onset_time./frame_time);  
CS_duration = info(5, 1);                           %in ms
CS_US_delay = info(6, 1);                           %in ms
US_onset_frame = round( (CS_onset_time + CS_US_delay)./frame_time );
US_duration = info(7, 1);                           %in ms

%loading behavioural response list and identifying animal as a learner or not
if set_type == 0        %trace
    mag_dir = [mag_an_direc 'trace\' mag_set_name];
elseif set_type == 1    %rand
    mag_dir = [mag_an_direc 'rand\' mag_set_name];
end

blink_list = load([mag_dir 'blink_list.txt']);
chance_blinks = load([mag_dir 'chance_blinks.txt']);
chance_blink_list = load([mag_dir 'chance_blink_vec.txt']);

bad_blinks = find(chance_blink_list > 3.5);
blink_list(bad_blinks) = 0;      %penalising for spontaneous blinking

mid_trial = round(no_trials./2);

if sum(blink_list) > 8
    learned = 1;
else
    learned = 0;
end