function [frame_time, zoom, n_chans, PMT_offsets, lpower] = SI_tif_info(stack_obj)
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

n_chansi = strfind(metadata, 'channelSave');
n_chans_newlinei = newlinesi - n_chansi;
n_chans_newlinei(n_chans_newlinei < 0) = nan;
[del, n_chans_newlinei] = nanmin(n_chans_newlinei);
n_chans_newlinei = newlinesi(n_chans_newlinei);
n_chans_string = metadata( (n_chansi + 14):(n_chans_newlinei - 1) );
eval(['n_chans = ' n_chans_string, ';']);
n_chans = length(n_chans);

%PMT channel offsets
PMT_offsetsi = strfind(metadata, 'channelOffsets');
PMT_offsets_newlinei = newlinesi - PMT_offsetsi;
PMT_offsets_newlinei(PMT_offsets_newlinei < 0) = nan;
[del, PMT_offsets_newlinei] = nanmin(PMT_offsets_newlinei);
PMT_offsets_newlinei = newlinesi(PMT_offsets_newlinei);
PMT_offsets = str2num(metadata( (PMT_offsetsi + 18):(PMT_offsets_newlinei - 2) ));

%checking if offsets were pre-subtracted
PMT_offsets_subi = strfind(metadata, 'SubtractOffsets');
PMT_offsets_sub_newlinei = newlinesi - PMT_offsets_subi;
PMT_offsets_sub_newlinei(PMT_offsets_sub_newlinei < 0) = nan;
[del, PMT_offsets_sub_newlinei] = nanmin(PMT_offsets_sub_newlinei);
PMT_offsets_sub_newlinei = newlinesi(PMT_offsets_sub_newlinei);
offset_str = metadata((PMT_offsets_subi + 19):(PMT_offsets_sub_newlinei - 2));
truei = findstr(offset_str, 'true');
%if offsets were pre-subtracted, forcing PMT_offsets to 0.
if isempty(truei) == 0
    PMT_offsets = [0, 0];
else
end

%getting laser power used
lpower_subi = strfind(metadata, 'hBeams.powers');
lpower_sub_newlinei = newlinesi - lpower_subi;
lpower_sub_newlinei(lpower_sub_newlinei < 0) = nan;
[del, lpower_sub_newlinei] = nanmin(lpower_sub_newlinei);
lpower_sub_newlinei = newlinesi(lpower_sub_newlinei);
lpower = metadata((lpower_subi + 16):(lpower_sub_newlinei - 1));
lpower = str2num(lpower);