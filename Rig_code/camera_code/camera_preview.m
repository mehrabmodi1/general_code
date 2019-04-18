clear all
close all

vid = videoinput('pointgrey');

figure('Toolbar','none',...
       'Menubar', 'none',...
       'NumberTitle','Off',...
       'Name','My Preview Window');

% Create the image object in which you want to display 
% the video preview data. Make the size of the image
% object match the dimensions of the video frames.

vidRes = vid.VideoResolution;
nBands = vid.NumberOfBands;
hImage = image( zeros(vidRes(2), vidRes(1), nBands) );

% Display the video data in your GUI.

preview(vid, hImage);

del = input('press enter when done with preview');

delete(vid)
close 'My Preview Window'