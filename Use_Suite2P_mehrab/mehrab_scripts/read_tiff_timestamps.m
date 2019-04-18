function [saved_times] = read_tiff_timestamps(direc)
%This function reads in the time stamps of each ScanImage tiff file in
%direc and saves them in a structure alongside the respective filenames.

%direc = 'D:\Data\Mehrab\20170906\fly2_axons\1\';

%making a list of tiff files in direc
old_dir = pwd;
cd(direc);
dir_contents = dir('*.tif');


%loop to read timestamp for each tif in direc
for tif_n = 1:size(dir_contents, 1)
    %reading in tiff meta data
    try
        curr_name = dir_contents(tif_n).name;
        img_obj = ScanImageTiffReader([direc, curr_name]);
    catch
        keyboard
    end
    
    meta_data = img_obj.descriptions;
    meta_data = meta_data{1};
    
    %----------------------------------------
    %reading in curr tiff's epoch and re-formatting the epoch string
    epochi = findstr(meta_data, 'epoch');
    t_string = meta_data( (epochi + 9):(epochi + 30) );                 %un-formatted time string saved by scanimage

    %re-formatting the string from the metadata to generate a matlab date object that can be used for time arithmetic.
    t_string_f = [t_string(1:4), '-', t_string(6:7), '-', t_string(9:10), ' '...
        t_string(12:13), ':', t_string(15:16), ':', t_string(18:19)];
    
    doti = findstr(t_string_f, '.');
    t_string_f(doti) = [];
    ep_time_obj = datetime(t_string_f, 'Format', 'yyyy-MM-dd hh:mm:ss');
    
    %----------------------------------------
    %reading in curr tiff's trigger time relative to epoch time and
    %creating an absolute time object for the current tiff.
    trigi = findstr(meta_data, 'TriggerTimestamps_sec');
    doti = findstr(meta_data(trigi:end), '.');
    trig_string = meta_data( (trigi + 24):(doti(1) - 2 + trigi) );
    trig_time_obj = ep_time_obj + seconds(str2num(trig_string));        %absolute time object for acq trigger for current tiff
    
    
    %----------------------------------------
    %keeping track of time-objs and filenames for all tifs in current
    %directory
    saved_times(tif_n).tstamp = trig_time_obj;
    saved_times(tif_n).name = curr_name;
    
end
