function [frame_time, zoom, n_chans] = SI_tif_info(stack_obj)
%syntax: [frame_time, zoom, n_chans] = SI_tif_info(stack_obj)

%%testing variables
% clear all 
% close all
% %stack_obj = ScanImageTiffReader('C:\Data\Data\Raw_data\20180409\cell2_OGB\buzz_trs_00001.tif');             %green alone
% stack_obj = ScanImageTiffReader('C:\Data\Data\Raw_data\20180419\cell1_nls_GC6f_seTau\brkin_00001.tif');     %green and red

metadata = stack_obj.metadata;
newlinesi = regexp(metadata, '[\n]');

zoomi = findstr(metadata, 'ZoomFactor');
zoom_newlinei = newlinesi - zoomi;
zoom_newlinei(zoom_newlinei < 0) = nan;
[del, zoom_newlinei] = nanmin(zoom_newlinei);
zoom_newlinei = newlinesi(zoom_newlinei);
zoom = str2num(metadata( (zoomi + 12):(zoom_newlinei - 1) ));


fr_ratei = strfind(metadata, 'scanFrameRate');
fr_rate_newlinei = newlinesi - fr_ratei;
fr_rate_newlinei(fr_rate_newlinei < 0) = nan;
[del, fr_rate_newlinei] = nanmin(fr_rate_newlinei);
fr_rate_newlinei = newlinesi(fr_rate_newlinei);
fr_rate = str2num(metadata( (fr_ratei + 15):(fr_rate_newlinei - 1) ));

avg_factori = strfind(metadata, 'logAverageFactor');
avg_factor_newlinei = newlinesi - avg_factori;
avg_factor_newlinei(avg_factor_newlinei < 0) = nan;
[del, avg_factor_newlinei] = nanmin(avg_factor_newlinei);
avg_factor_newlinei = newlinesi(avg_factor_newlinei);
avg_factor = str2num(metadata( (avg_factori + 18):(avg_factor_newlinei - 1) ));

fr_rate = fr_rate./avg_factor;
frame_time = 1./fr_rate;      %in s
frame_time = frame_time.*1000; % in ms


n_chansi = strfind(metadata, 'channelSave');
n_chans_newlinei = newlinesi - n_chansi;
n_chans_newlinei(n_chans_newlinei < 0) = nan;
[del, n_chans_newlinei] = nanmin(n_chans_newlinei);
n_chans_newlinei = newlinesi(n_chans_newlinei);
n_chans_string = metadata( (n_chansi + 14):(n_chans_newlinei - 1) );
eval(['n_chans = ' n_chans_string, ';']);
n_chans = length(n_chans);