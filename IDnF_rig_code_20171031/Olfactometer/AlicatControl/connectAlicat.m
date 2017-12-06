function varargout=connectAlicat(aliComm)
% function aliComm=connectAlicat(aliComm)
%
% Form a connection to the Alicat controller. If an argument is given then
% the connection to the specified device is closed.
%
%
% Rob Campbell - 20th March 2008 - CSHL

com_port = 'COM12';

if nargin==0
    try
        global aliComm
        aliComm=serial(com_port,...
            'TimeOut', 2,...
            'BaudRate', 19200,...
            'Terminator','CR');

        fopen(aliComm)
        set(aliComm,'Tag','Alicat_serial');
        varargout{1}=aliComm;
    catch
        obj_n = 1;
        %making sure that MFC serial port is closed to be opened afresh.
        while obj_n ~= 0
            a = instrfind;
            disp(a)
            obj_n = input('Which item in the list is the MFC serial port? Enter 0 if absent.');
            
            if obj_n > 0
                fclose(a(obj_n))
                delete(a(obj_n))
            else
            end
        end
        
        global aliComm
        aliComm=serial(com_port,...
            'TimeOut', 2,...
            'BaudRate', 19200,...
            'Terminator','CR');

        fopen(aliComm)
        set(aliComm,'Tag','Alicat_serial');
        varargout{1}=aliComm;
        
    end
    
elseif nargin==1
    fclose(aliComm)
    delete(aliComm)
    clear global aliComm
end
