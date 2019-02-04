% Example script to run a simple pRF analysis on a voxel-by-voxel basis
% and view the results.
%
% Neuronal model:      Single Gaussian function
% Receptive field:     Circular (isotropic)
% Input specification: Polar coordinates


%% Settings
clc; close all;

% Check whether parent directory was provided when calling MATLAB from
% command line. Check if it was not provided. 
try
    par_dir;
catch
    par_dir = '/media/sf_D_DRIVE/MotDepPrf';
end

% Root directory for glm and data
root_dir = fullfile(par_dir, 'Analysis','S02','06_bayesPrf');

% Directory with ROI definitions
surf_dir = fullfile(par_dir, 'Analysis','S02','Anatomy', 'Session02');

% Directory for creating GLM
glm_dir = fullfile(root_dir,'pRF_results');

% Set VOI name
try
  voi_name;
catch
  voi_name = 'S02_2H_allRois_1mm';
end

% Set model name
try
  mdl_name;
catch
  mdl_name = 'spm_prf_fcn_gaussian_polar';
end

% Motion condition
try
  mtn_cnd;
catch
  mtn_cnd = '_motLoc';
end

% Repetition time
TR = 2;
% Echo time
TE = 0.020;
% Bins per TR (default 16 or number of slices)
nmicrotime    = 35;
% Duration of stimuli (secs)
stim_duration = 1.5;
% Diameter of stimuli in degrees
stim_diameter = 24;

%% Derive and adjust settings

% Update glm_dir
glm_dir = fullfile(glm_dir, [voi_name mtn_cnd]);

% Which sessions to include
sess = 1:4;
num_sess = length(sess);

%% Prepare inputs

% Load exp info
load(fullfile(root_dir,'expInfo',['ApnFrm',mtn_cnd,'.mat']));
U = prepare_inputs_polar_samsrf(ApFrm,TR, nmicrotime, stim_duration, stim_diameter);
% Remove empty fields from structure
empty_elems = arrayfun(@(s) all(structfun(@isempty,s)), U);
U(empty_elems) = [];

% The timeseries from each session are stored in a VOI_xx.mat file. Build a
% cell array of the VOI files for each session.
xY = cell(1,num_sess);
for i = 1:num_sess
    filename = sprintf('VOI_%s_%d.mat',voi_name, i);
    xY{i}    = fullfile(glm_dir,filename);
end

%% Specify pRF model

% Load SPM for timing information / image dimensions
SPM = load(fullfile(glm_dir,'SPM.mat'));
SPM = SPM.SPM;

% Update SPM path as we don't know where this example will be saved
SPM.swd = glm_dir;

% Set pRF specification options
options = struct('TE',TE, ...             % echo time
                 'voxel_wise',true, ...   % per voxel (true) or ROI (false)
                 'model', mdl_name,... % pRF function (spm_prf_fcn...)
                 'hE',6, ...              % expected log precision of noise
                 'P',[], ...              % starting parameters
                 'B0',7, ...              % fMRI field strength (teslas)
                 'avg_sess',true, ...     % average sessions' timeseries
                 'avg_method','mean',...  % accepts 'mean' or 'eigen'
                 'name', [voi_name mtn_cnd]);


% Specify pRF model (.mat file will be stored in the GLM directory)
PRF = spm_prf_analyse('specify',SPM,xY,U,options);
num_voxels = size(PRF.Y.y, 2);

%% Estimate all voxels (slow)
voxel = 1:num_voxels;

% Model to estimate
prf_file = fullfile(glm_dir,['PRF_' voi_name mtn_cnd '.mat']);

% Estimation options
options = struct('init', 'GLM_P', ...    % Initialization.
                 'use_parfor', 6, ...% Parallelization
                 'nograph', true, ...    % If true, disables plots
                 'voxels', voxel ...     % Voxels indices (optional)
                 );

% Estimate
PRF_est = spm_prf_analyse('estimate',prf_file,options);

% Review
spm_prf_review(prf_file);
