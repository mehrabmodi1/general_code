function FlipValve_EP(lines,state)
%Flips valves through USB6501
%
%Input
%   -lines: determines which line to flip.  You can determine with number
%   (Index number of USB6501 lines) or character of valve names (i.e. 'vial1', 'vial2', 'NO')
%   Or you can use 'all' to flip all valves
%   You can also define multiple valves with cell array.
%
%   -state: 0 or 1.  Caution: 0 is energized state
%


%Only connect to the board the first time the function is called.  
persistent vial_switch olfSS ValveState; %see connectToUSB6501_EP.m

if isempty(olfSS)
    [vial_switch,olfSS,ValveState]=connectToUSB6501_EP;
end

if isnumeric(lines)
    ValveState(lines)=state;
elseif ischar(lines) && strcmpi(lines,'all')
    ValveState(1:end)=state;
elseif ischar(lines)
    lines={lines};
end
if iscell(lines)
    if isscalar(state)
        state=state*ones(size(lines));
    end
    for i=1:length(lines)
        if strcmpi(lines{i},'vial9')
            lines{i}='NO';
        end
        ValveState(strcmpi(lines{i},vial_switch(:,1)))=state(i);
    end
end
outputSingleScan(olfSS,ValveState);