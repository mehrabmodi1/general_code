%function [frame_size, n_frames, im_mat] = load_binary_im(direc, file_n)
%this function makes a list of the loadable .bin image files and then loads
%the file numbered file_n in this list into memory.
direc = 'D:\Data\CSHL\Resonant\20161004_2\zstack\zoomed\';
file_n = 1;

orig_dir = pwd;
cd(direc);
bin_list = dir('*.bin');
n_files = size(bin_list, 1);

if file_n > n_files 
    error('file number specified greater than the number of bin files in directory.')
else
end

f_name = bin_list(file_n).name;
[framesize,frames] = Binary2Image(f_name);

cd(orig_dir);
%end