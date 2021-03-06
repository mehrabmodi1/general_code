function [time_stamp] = parse_tiff_timestamp(img_obj)
%syntax: [time_stamp] = parse_tiff_timestamp(stack_obj)
%This function parses the image description read out by ScanImageTiffReader
%to give the time-stamp used later in analysis to align stim params to
%acquired imaging data.


meta_data = img_obj.descriptions;
meta_data = meta_data{1};

%reading in curr tiff's epoch and re-formatting the epoch string
epochi = findstr(meta_data, 'epoch');
t_string = meta_data( (epochi + 9):(epochi + 30) );                 %un-formatted time string saved by scanimage

%re-formatting the string from the metadata to generate a matlab date object that can be used for time arithmetic.
t_string_f = [t_string(1:4), '-', t_string(6:7), '-', t_string(9:10), ' '...
    t_string(12:13), ':', t_string(15:16), ':', t_string(18:19)];

doti = findstr(t_string_f, '.');
t_string_f(doti) = [];
ep_time_obj = datetime(t_string_f, 'Format', 'yyyy-MM-dd hh:mm:ss');

%reading in curr tiff's trigger time relative to epoch time and
%creating an absolute time object for the current tiff.
trigi = findstr(meta_data, 'TriggerTimestamps_sec');
doti = findstr(meta_data(trigi:end), '.');
trig_string = meta_data( (trigi + 24):(doti(1) - 2 + trigi) );
trig_time_obj = ep_time_obj + seconds(str2num(trig_string));        %absolute time object for acq trigger for current tiff
time_stamp = trig_time_obj;