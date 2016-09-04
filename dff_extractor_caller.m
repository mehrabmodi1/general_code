%This script calls dff_extractor, reading dataset paths from a list of
%paths in a text file, so as to extract data from multiple datasets
%automatically.

list_path = 'C:\Data\CSHL\dataset_list_20141208.txt';

fid=fopen(list_path);
counter = 0;
while 1
    counter = counter + 1;
    direc = fgetl(fid);
    if ~ischar(direc),   break,   end
    direc = [direc '\'];
    
    
    disp(['dataset number ' int2str(counter)])
    dff_extractor(direc, 0)
    
end
fclose(fid);