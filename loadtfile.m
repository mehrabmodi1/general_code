function [t] = loadtfile (fname)
%loads t-files written by MClust, containing spike time-stamps. the input
%fname is the path to the t-file, including its filename.

    fid = fopen(fname , 'r', 'b');
    
    hdr = fread(fid, 400, '*char')';
    fseek (fid, findstr(hdr, '%%ENDHEADER') + 11, 'bof');
    t = fread(fid, inf, 'uint32');
    fclose (fid);