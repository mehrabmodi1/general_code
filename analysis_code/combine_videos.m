clear all
close all

path_vec = 'C:\Data\Data\Analysed_data\Analysis_results\air_gap_traj_an\vert_aligned\';

t_subpaths = [{'0s'}, {'25s'}];
od_subpaths = [{'prd'}, {'unprd'}];
vid_counter = 0;
%reading in all video frames
for time_n = 1:2
    time_i = t_subpaths{time_n};
    for od_n = 1:2
        vid_counter = vid_counter + 1;
        od_i = od_subpaths{od_n};
        curr_path = [path_vec, time_i, ' ', od_i, '.mj2'];
        eval(['v', num2str(vid_counter), '= VideoReader(curr_path)']); 
        keyboard
    end
end

while hasFrame(v1)
    frame1 = readFrame(v1);
    frame2 = readFrame(v2);
    frame3 = readFrame(v3);
    frame4 = readFrame(v4);
    keyboard
end
