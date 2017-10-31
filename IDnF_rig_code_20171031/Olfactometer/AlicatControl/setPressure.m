function setPressure(AC,MFC,pressure)
% function setPressure(AC,MFC,pressure)
%
% Set a given MFC to maintain a fixed pressure. 
%
% Inputs
% AC - serial port object
% MFC - string specifying the controller ID.
% pressure - desired pressure in PSI
%
% Rob Campbell - June 2010


writeLoop(AC,MFC,3)

units=(pressure/160)*64E3;

fprintf(AC,[MFC,num2str(round(units))])
fscanf(AC);