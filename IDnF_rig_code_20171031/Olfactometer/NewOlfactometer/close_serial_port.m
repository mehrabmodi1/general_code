function [] = close_serial_port(port_num)
%This function uses instrfind to check if the serial port named in
%port_string is open already and closes it so that it can be opened
%gracefully if necessary.

port_list = instrfind;

if isnumeric(port_num) == 1
    port_num = num2str(port_num);
else
end

port_string = ['Serial-COM', port_num];

port_n = [];
%identifying relevant port object in port_list
for open_port_n = 1:size(port_list, 2)
    curr_name = port_list(open_port_n).Name;
    if strcmp(curr_name, port_string) == 1
        port_n = [port_n, open_port_n];     %in case of multiple instantiations
    else
    end
end

%closing all instantiations of specified port
if length(port_n) >= 1
    for close_port_n = 1:length(port_n)
        curr_port_n = port_n(close_port_n);
        fclose(port_list(curr_port_n));
        delete(port_list(curr_port_n));
        clear port_list(curr_port_n)
    end
else
end