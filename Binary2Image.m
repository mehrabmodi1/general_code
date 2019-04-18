%% load one or more frames from binary image file saved on photoncerber
%% call the function from the folder that contains the .bin file

function [framesize,frames,ImageOut] = Binary2Image(datafile,numframes);
    [~,datafilename] = fileparts(datafile);
    byteorder = 'ieee-le';
    prec = 'int16';
    bytesize = 2;
    
    % open file
    fid = fopen(datafile,'r',byteorder);
    fseek(fid,0,'bof'); % Reset the cursor
    
    % read frame size
    framesize = fread(fid,2,prec,0,byteorder);
    pixels = framesize(1)*framesize(2);
        
    % get total number of frames
    s = dir(datafile);
    frames = (s.bytes - 4)/(pixels*bytesize); % 2 bytes per value
    
    if nargin==2
        frames = numframes;
    end
    
    % read frame 1 and create the tiff stack
    imwrite(uint16(fread(fid,[framesize(2) framesize(1)],prec,0,byteorder)),...
        [datafilename,'.tif']);
    
    % append to the stack
    for f = 2:frames
        %ImageOut(:,:,f) = uint16(fread(fid,[framesize(2) framesize(1)],prec,0,byteorder));
        imwrite(uint16(fread(fid,[framesize(2) framesize(1)],prec,0,byteorder)),...
        [datafilename,'.tif'],'WriteMode','append');
    end
    
    fclose(fid);
end

    
    