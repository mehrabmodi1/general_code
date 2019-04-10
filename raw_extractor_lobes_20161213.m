clear all
close all


direc_list = 'D:\Data\CSHL\dataset_list_sustained_OK107xsytGCaMP6s_20161212.txt';

fid = fopen(direc_list);
direc_counter = 0;

while 1
        
    direc = fgetl(fid);

    if ischar(direc) ~= 1
        break
    else
    end

    direc_counter = direc_counter + 1;

    keyboard 
        
        
        
end