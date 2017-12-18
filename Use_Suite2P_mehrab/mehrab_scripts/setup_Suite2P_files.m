function [ft_factor] = setup_Suite2P_files(raw_direc_base, raw_direc, m_db_file_direc)
%This function copies over the master_file.m and make_db.m files needed to
%run Suite2P to the specified [raw_direc_base raw_direc] from the specified m_db_file_direc.
%It also edits both files according to the dataset in the current
%raw_direc.

    %creating make_db and master_file as needed to run Suite2P
    prev_direc = pwd;
    cd([raw_direc_base, raw_direc])
    bslashi = findstr(raw_direc, '\');
    
    %This function creates a make_db.m file from scratch for each dataset.
    %This typically needs to be edited for every dataset.
    fid = fopen('make_db.m', 'w');
    fprintf(fid, ['i = 1;', '\n']);
    fprintf(fid, ['db(i).mouse_name    = ''' raw_direc(1:bslashi(1) ) ''';' '\n']);
    fprintf(fid, ['db(i).date          = ''' raw_direc((bslashi(1)+1):bslashi(2)) ''';' '\n']);
    fprintf(fid, ['db(i).expts         = [''' raw_direc((bslashi(2)+1):bslashi(3)) '''];' '\n']);
    fprintf(fid, 'db(i).nchannels     = 1;\n');
    fclose(fid);
    
    %This funtion copies over master_file.m from the path given to it as
    %one of it's inputs. This typically needs to be edited very rarely.
    copyfile(m_db_file_direc, [raw_direc_base, raw_direc]);
    
    %reading in frame rate and editing local master_file to reflect it
    dir_contents = dir('*.tif');
    dir_contents = dir_contents.name;
    reader = ScanImageTiffReader(dir_contents);
    desc = reader.metadata;
    base_ft_i = findstr(desc, 'SI.hRoiManager.scanFramePeriod');
    base_frame_time = desc(1, (base_ft_i + 33):(base_ft_i + 40));       %in s, frame time before online averaging
    base_frame_time = str2num(base_frame_time);
    ft_factori = findstr(desc, 'SI.hScan2D.logAverageFactor');
    newlinei = findstr(desc(ft_factori:end), sprintf('\n'));
    ft_factor = desc(1, (ft_factori + 30):(ft_factori + newlinei(1) - 1));
    ft_factor = str2num(ft_factor);
    frame_time = base_frame_time.*ft_factor;
    
    %editing masterfile with parameters for current dataset
    master_file = regexp(fileread('master_file.m'), '\n', 'split');
    fr_line = master_file{38};
    
    frame_time_factor = num2str(1./frame_time, 2);
    try
        fr_line(1, 31:(31 + length(frame_time_factor) - 1) ) = frame_time_factor;
    catch
        keyboard
    end
    master_file{38} = fr_line;
    
    fid = fopen('master_file.m', 'w');
    fprintf(fid, '%s\n', master_file{:});
    fclose(fid);
    
    cd(prev_direc)


end