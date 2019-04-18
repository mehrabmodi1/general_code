 function [] = z_sligner_res()
% 
% hSI = evalin('base','hSI');        % get hSI from the base workspace
% 
% %trial loop
% while 1
%     
%     %check SI status loop
%     SI_status = 0;
%     while SI_status == 0
%         curr_state = hSI.acqState;
%         SI_status = strcmp(curr_state, 'loop_wait');
%         pause(5)        
%     end
%     
    
    
    %loading z_stack
    z_stack_mat = load_z_stack();
    
    %loading last saved trial tif
    aq_direc = curr_aq_direc;
    dir_contents = dir([aq_direc, '*.tif']);
    n_files = size(dir_contents, 1);
   
    
    %identifying plane for tr 1
    if n_files == 1
        tif_stack = load_tif_stack([aq_direc, dir_contents(1).name]);
        frame_ave = mean(tif_stack, 3);
        
        plane_matches = compare_with_zstack(frame_ave, z_stack_mat, 0, 1);
        [del, plane_n] = max(plane_matches);
        
        %PICK UP THREAD HERE:
        %currently failing to accuratly identify the matching plane in the
        %z-stack. Unclear why correlation is not a good measure to compare
        %two images.
       
        
    else
    end
    
    
   
    %identifying latest tif file
    max_date = [dir_contents(1).datenum, 1];
    if n_files > 1
        for file_n = 2:n_files
            curr_date = dir_contents(file_n).datenum;
            if curr_date > max_date(1, 1);
                max_date = [curr_date, file_n];
                
            else
            end
        end
    else
    end
    
    keyboard
    
%end




end