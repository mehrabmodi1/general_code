function [direc] = dataset_assembler(direc, dir_counter)
%dataset_assembler syntax: direci = dataset_assembler(direc)
%dataset_assembler looks for other folders at the same level as  direc,
%of the form [direc '1.x'] where x = 1, 2, 3..., beginning with 1. If such
%folders are found, all .tif and .xls files in direc, [direc '1.1'], [direc
%'1.2']... will be copied into a new folder called [direc '_assembled']


%D:\Data\data\BlinkView\20120113\field2_rand_final_final-ROI

%loop to look for fragment directories
i = 0;
count = 0;

while i == 0
    count = count + 1;
    checkdir = [direc(1, 1:(end - 4)) '1.' int2str(count) '-ROI'];
    checkdiri = [direc(1, 1:(end - 4)) '_assembled\'];
    
    if isdir(checkdir) == 1
        
        if isdir(checkdiri) == 1        %checking for existence of _assembled directory and then quitting dataset_assembler
            i = 0;
            count = 0;
            direc = (checkdiri(1, 1:(end-1)));
            
        else

            i = 0;
            %disp('yes')
        end
        
    else
        i = 1;
        count = count - 1;
        %disp('no')

    end
end

if count > 0
    dir_assembled = [direc(1, 1:(end - 4)) '_assembled\' ];
    mkdir(dir_assembled);
    
    %loop to open fragment directories
    trial_no = 0;
    rand_timings_listi = [];
    for dirno = 0:count
        if dirno == 0
            direci = direc;
            
        else
            direci = [direc(1, 1:(end - 4)) '1.' int2str(dirno) '-ROI'];
        end
        

        no_trials = floor(( length(dir(direci)) - 2 )./2);

        %loading RandCStimings if any
        dir_cont_xls = dir([direci '\*.xls']);
        dir_cont_tif = dir([direci '\*.tif']);
        
        if (length(dir_cont_xls) - length(dir_cont_tif) ) == 1
            rand_timings_list = load([direci '\rand CS timings.xls']); 
            rand_timings_list((no_trials+1):end) = [];
            rand_timings_listi = [rand_timings_listi, rand_timings_list];
        else
        end


        for curr_trial_no = 1:no_trials
            %copying .tiff and .xls files
            copyfile([direci '\Trial - ' int2str(curr_trial_no-1) '-1.tif'], [dir_assembled 'Trial - ' int2str(trial_no) '-1.tif' ]);
            copyfile([direci '\Trial - ' int2str(curr_trial_no-1) '-1.xls'], [dir_assembled 'Trial - ' int2str(trial_no) '-1.xls' ]);
            trial_no = trial_no + 1;
            
%             if dir_counter == 11
%                 keyboard
%             else
%             end
            
                       
        end
       
        
        
        
    end
    save([dir_assembled 'rand CS timings.xls'],'rand_timings_listi','-ASCII');
    
    direc = dir_assembled;
else
    direc = [direc '\'];
end
