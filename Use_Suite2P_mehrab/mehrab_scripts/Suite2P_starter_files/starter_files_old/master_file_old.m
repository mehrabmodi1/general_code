%% SET ALL DEFAULT OPTIONS HERE
% check out the README file for detailed instructions (and extra options)
%addpath('D:\CODE\MariusBox\runSuite2P') % add the path to your make_db file

% overwrite any of these default options in your make_db file for individual experiments
make_db; % RUN YOUR OWN MAKE_DB SCRIPT TO RUN HERE

ops0.toolbox_path = 'D:\Data\Bitbucket_repos\Suite2P_copy\';
if exist(ops0.toolbox_path, 'dir')
	addpath(genpath(ops0.toolbox_path)) % add local path to the toolbox
else
	error('toolbox_path does not exist, please change toolbox_path');
end

% mex -largeArrayDims SpikeDetection/deconvL0.c (or .cpp) % MAKE SURE YOU COMPILE THIS FIRST FOR DECONVOLUTION

ops0.useGPU                 = 1; % if you can use an Nvidia GPU in matlab this accelerates registration approx 3 times. You only need the Nvidia drivers installed (not CUDA).

% root paths for files and temporary storage (ideally an SSD drive. my SSD is C:/)
ops0.RootStorage            = 'D:\Data\CSHL\Resonant\'; % Suite2P assumes a folder structure, check out README file
ops0.temp_tiff              = 'D:\Data\Suite2P_scratch\temp.tif'; % copy each remote tiff locally first, into this file
ops0.RegFileRoot            = 'D:\Data\Suite2P_scratch\';  % location for binary file
ops0.DeleteBin              = 1; % set to 1 for batch processing on a limited hard drive
ops0.ResultsSavePath        = 'D:\Data\Suite2P_results\'; % a folder structure is created inside
ops0.RegFileTiffLocation    = ['D:\Data\Suite2P_registered\']; %'D:/DATA/'; % leave empty to NOT save registered tiffs (slow)

% registration options
ops0.doRegistration         = 1; % skip (0) if data is already registered
ops0.showTargetRegistration = 1; % shows the image targets for all planes to be registered
ops0.PhaseCorrelation       = 1; % set to 0 for non-whitened cross-correlation
ops0.SubPixel               = Inf; % 2 is alignment by 0.5 pixel, Inf is the exact number from phase correlation
ops0.NimgFirstRegistration  = 500; % number of images to include in the first registration pass 

% cell detection options
ops0.clustModel             = 'neuropil'; % standard or neuropil
ops0.neuropilSub            = 'model'; % none, surround or model
ops0.ShowCellMap            = 1; % during optimization, show a figure of the clusters
ops0.Nk0                    = 500; % how many clusters to start with
ops0.Nk                     = 200;  % how many clusters to end with (before anatomical segmentation)
ops0.sig                    = 0.5;  % spatial smoothing length in pixels; encourages localized clusters
ops0.nSVDforROI             = 250; % how many SVD components for cell clustering
ops0.NavgFramesSVD          = 5000; % how many (binned) timepoints to do the SVD based on
clustrules.diameter         = 20; % expected diameter of cells (used for 0.25 * pi/4*diam^2 < npixels < 10*pi/4*diam^2)

% spike deconvolution options
ops0.imageRate              = 50;   % imaging rate (cumulative over planes!). Approximate, for initialization of deconvolution kernel.
ops0.sensorTau              = 2; % decay half-life (or timescale). Approximate, for initialization of deconvolution kernel.
ops0.maxNeurop              = Inf; % for the neuropil contamination to be less than this (sometimes good, i.e. for interneurons)
ops0.recomputeKernel        = 1; % whether to re-estimate kernel during optimization (default kernel is "reasonable", if you give good timescales)
ops0.sameKernel             = 1; % whether the same kernel should be estimated for all neurons (robust, only set to 0 if SNR is high and recordings are long)

db0 = db;
%% RUN THE PIPELINE HERE
for iexp = 1 %:length(db)
    run_pipeline(db(iexp), ops0, clustrules);
    
    % deconvolved data into (dat.)cl.dcell, and neuropil subtraction coef
    add_deconvolution(ops0, db0(iexp), clustrules);
    
    % add red channel information (if it exists)
    %     run_REDaddon(iexp, db, ops0) ;
end
%% STRUCTURE OF RESULTS FILE
% 
% cell traces are in dat.F.Fcell
% neuropil traces are in dat.F.FcellNeu
% neuropil subtraction coefficient is dat.cl.dcell{i}.B(3)
% baseline is dat.cl.dcell{i}.B(2)
% anatomical cell criterion is in dat.cl.isroi
% manual overwritten cell labels are in dat.cl.iscell
% dat.cl.dcell{i}.st are the deconvolved spike times (in frames)
% dat.cl.dcell{i}.c  are the deconvolved amplitudes
% dat.cl.dcell{i}.kernel is the estimated kernel