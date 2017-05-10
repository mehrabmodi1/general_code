results_direc = 'D:\Data\Suite2P_results\';
reg_tif_direc = 'D:\Data\Suite2P_registered\';
raw_direc = 'D:\Data\CSHL\Resonant\20170313\MB149BX20xsytGC6f - Copy\odor_trials\';

%% Pre-prep
%Getting rid of empty .tiff files in raw data folder(s) left behind by ScanImage
remove_small_tifs(raw_direc);


%% Running GUI_Suite2P_Main
gui_h = GUI_Suite2P_Main;
disp('close GUI window once its done to continue')
uiwait(gui_h)       %wait until GUI window is closed

curr_set = 'Suite2Ptest\20160517\1\';


%% Reading in registered tiffs and displaying for manual removal of trials with z-drift
prev_direc = pwd;
cd([ref_tif_direc, curr_set])

keyboard
cd(prev_direc)
%% Loading in Suite2P analysis results file
prev_direc = pwd;
cd([results_direc, curr_set])
dir_contents = dir;
dir_contents(1:2) = [];
n_files = size(dir_contents, 1);
%finding most recent file
max_datenum = [];
for file_n = 1:n_files
    fname = dir_contents(file_n).name;
    if isempty(findstr(fname, 'Nk200')) == 1
        continue
    else
    end
    curr_datenum = dir_contents(file_n).datenum;
    if isempty(max_datenum) == 1
        max_datenum = [curr_datenum, file_n];
    elseif isempty(max_datenum) == 0
        if max_datenum(1) < curr_datenum == 1
            max_datenum(1) = curr_datenum;
            max_datenum(2) = file_n;
        else
        end
    end
    
end

data_mat = load([results_direc, curr_set, dir_contents(max_datenum(2)).name]);
data_mat = data_mat.dat;

%% Creating ROI matrix that includes only ROIs classified as real cells
list_rois = data_mat.cl.isroi;
list_cell_rois = find(list_rois == 1); 
n_rois = length(list_cell_rois);
size_im = size(data_mat.cl.k1);
ROI_mat = zeros(size_im(1), size_im(2), n_rois);        %each ROI separated along third dimension.
for roi_n = 1:n_rois
    roi_n_orig = list_cell_rois(roi_n);
    roi_pix = data_mat.stat(roi_n_orig).ipix;           %pixels for current roi
    ROI_im = zeros(size_im(1), size_im(2));
    ROI_im(roi_pix) = 1;
    ROI_mat(:, :, roi_n) = ROI_im;
end

%% Extracting raw fluorescence traces from registered tiff images
cd([reg_tif_direc, curr_set, 'Plane1\']);
dir_contents = dir;
%PICK UP THREAD HERE
stack=ScanImageTiffReader([raw_direc, 'KC_prep_00028.tif']).data();
m_data=ScanImageTiffReader([raw_direc, 'KC_prep_00028.tif']);
m_data = m_data.metadata;
keyboard



%- read in tiff files in order of trials, extract raw traces for each roi,
%write extracted traces to disk
%- downsample registered tiff files by averaging frames and overwrite
%original files to save disk space

