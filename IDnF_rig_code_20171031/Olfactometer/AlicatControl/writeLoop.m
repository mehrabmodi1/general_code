function varargout=writeLoop(AC,MFC,loop)
% function loop=writeLoop(AC,MFC,loop)
%
% Set the loop mode of an MFC
% 
% Inputs
% AC - serial port object 
% MFC - string specifying the controller ID
% loop - a scalar. 1 for mass, 2 for volume, 3 for pressure. 
%
% Outputs
% loop - [optional] string confirming the new loop type. 
%
% Rob Campbell June 2010

setFlow(AC,0,MFC) %Make sure we're zero before doing anything

%9460 is mass flow (or should be!)
L(1)=9460;
L(2)=L(1)-2^10+768;
L(3)=L(1)-2^10+256;

fprintf(AC,[MFC,'$$W20=',num2str(L(loop))]);
fscanf(AC);

if nargout==1
    varargout{1}=readLoop(AC,MFC);    
end
