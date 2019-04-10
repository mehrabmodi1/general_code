function varargout=readLoop(AC,MFC)
% function loop=readLoop(AC,MFC)
%
% Query the loop mode of an MFC
% 
% Inputs
% AC - serial port object 
% MFC - string specifying the controller ID
%
% Outputs
% loop - string giving the loop type. If no output arg then 
% just print to screen. 
%
% Rob Campbell June 2010


fprintf(AC,[MFC,'$$R20']);
F=fscanf(AC);


loop=regexp(F,'= (.*)','tokens');
loop=str2num(loop{1}{1});

%9460 is mass flow (or should be!)
mass=9460;
volumetric=mass-2^10+768;
pressure=mass-2^10+256;

switch loop
    case volumetric
        loop='volumetric';
    case mass
        loop='mass';
    case pressure
        loop='pressure';
end


if nargout==0
    fprintf('%s\n',loop)
end
if nargout>0
    varargout{1}=loop;
end
