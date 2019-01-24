%% Runs GLM analysis and extracts timeseries

% Settings
clc;close all;
spm('defaults','FMRI')

% Directory into which to download example dataset
data_root_dir = '/media/sf_D_DRIVE/MotDepPrf/Analysis/S02/06_bayesPrf/';
data_dir      = '/media/sf_D_DRIVE/MotDepPrf/Analysis/S02/06_bayesPrf/data/';
surf_dir      = '/media/sf_D_DRIVE/MotDepPrf/Analysis/S02/Anatomy/Session02/';

% Directory for creating GLM
glm_dir  = '/media/sf_D_DRIVE/MotDepPrf/Analysis/S02/06_bayesPrf/pRF_results/';

% Set VOI name
try
  voi_name;
catch
  voi_name = 'S02_2H_V1_1mm';
end

% Motion condition
try
  mtn_cnd;
catch
  mtn_cnd = '_expn';
end

% Stepping direction
try
  stp_drc;
catch
  stp_drc = '_outward';
end

% Repetition time
TR = 2;
% Bins per TR (default 16 or number of slices)
nmicrotime    = 35;
% Duration of stimuli (secs)
stim_duration = 1.4;
% Diameter of stimuli in degrees
stim_diameter = 17;

%% Derive and adjust settings

% Change the current folder to the folder of this m-file.
if(~isdeployed)
    tmp = matlab.desktop.editor.getActive;
    cd(fileparts(tmp.Filename));
end

% Set present directory to variable
start_dir = pwd;

% Update glm_dir
glm_dir = fullfile(glm_dir, [voi_name mtn_cnd stp_drc]);

% Derive vector for inward or outward runs
switch stp_drc
    case '_outward'
        sess = 1:2:11;
    case '_inward'
        sess = 2:2:11;
end

% Derive number of sessions
nsess = length(sess);

%% Prepare onsets

% load exp info
load(fullfile(data_root_dir,'expInfo',['ApnFrm',mtn_cnd,stp_drc,'.mat']));
U = prepare_inputs_polar_samsrf(ApFrm,TR, nmicrotime, stim_duration, stim_diameter);
% remove empty fields from structure
empty_elems = arrayfun(@(s) all(structfun(@isempty,s)), U);
U(empty_elems) = [];

bins_d = linspace(3.4, 8.5, 8);

% Build time x screen bins matrix (3x3 screen bins)
onsets_matrix = zeros(length(U), length(bins_d));
for t = 1:length(U)
    % Loop over pixels activated at this time point
    for activated_pixel = 1:length(U(t).dist)
        % Get location
        dist  = U(t).dist(activated_pixel);
        angle = U(t).angle(activated_pixel);

        % Identify closest bin
        [~,bin_idx] = min(abs(bins_d-dist));

        onsets_matrix(t,bin_idx) = onsets_matrix(t,bin_idx) + 1;
    end
end

% Remove empty bins
onsets_matrix = onsets_matrix(:,any(onsets_matrix));
num_regressors = size(onsets_matrix,2);

% SPM inputs
names = cell(1,num_regressors);
onsets = cell(1,num_regressors);
durations = cell(1,num_regressors);

for t = 1:num_regressors
    names{t} = ['Bin' num2str(t)];
    onsets{t} = (find( onsets_matrix(:,t) ) - 1) * TR;
    durations{t} = 0;
end

% Make output directory
if ~exist(glm_dir,'file')
    mkdir(glm_dir);
end

save(fullfile(glm_dir, ['onsets_',voi_name, mtn_cnd,stp_drc,'.mat']), 'names', 'onsets', 'durations');

%% Specify first level design

% Load generic matlabbatch for fmri_spec, fmri_est and con(trast)
load('first_level_batch.mat');

% Session-specific options
for i = 1:nsess
    str_run_num = sprintf('%02d', sess(i));
    str_run = ['func', str_run_num, '_SlTiSPM_MoCoSPM_msk', mtn_cnd, '.nii'];
    epis     = spm_select('ExtFPList',data_dir, str_run, 1:999);

    matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans     = cellstr(epis);
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi     = cellstr(fullfile(glm_dir, ['onsets_',voi_name, mtn_cnd,stp_drc,'.mat']));
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf = 112;
end

% Model spec options
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = nmicrotime;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(glm_dir);
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'none';

% Specify spm.mat location
matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(fullfile(glm_dir,'SPM.mat'));
matlabbatch{3}.spm.stats.con.spmmat = cellstr(fullfile(glm_dir,'SPM.mat'));

% Initialise job configuration
spm_jobman('initcfg')
% Run job
spm_jobman('run',matlabbatch);

cd(start_dir);

%% Build a mask of voxels which survive p < 0.001
clear matlabbatch;
matlabbatch{1}.spm.stats.results.spmmat = cellstr(fullfile(glm_dir,'SPM.mat'));
matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'none';
matlabbatch{1}.spm.stats.results.conspec.thresh = 0.005;
matlabbatch{1}.spm.stats.results.conspec.extent = 0;
matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
matlabbatch{1}.spm.stats.results.units = 1;
matlabbatch{1}.spm.stats.results.print = false;
matlabbatch{1}.spm.stats.results.write.tspm.basename = 'mask_uncorrected';

% Run job
spm_jobman('run',matlabbatch);

cd(start_dir);

%% Extract timeseries from surface voxels which survive p < 0.001

% Identify masks
mask1   = fullfile(surf_dir,[voi_name, '.nii.gz']);
[p1,p2,~] = fileparts(mask1);
mask2   = fullfile(glm_dir,'spmF_0001_mask_uncorrected.nii');

% gunzip
unix(['gunzip ', mask1]);

% Prepare batch
load('extract_timeseries_batch.mat');

matlabbatch{1}.spm.util.voi.name = [p2(1:end-4)];
matlabbatch{1}.spm.util.voi.spmmat = cellstr(fullfile(glm_dir,'SPM.mat'));
matlabbatch{1}.spm.util.voi.adjust = 1;
matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold = 0.5;
matlabbatch{1}.spm.util.voi.roi{1}.mask.image = cellstr(mask1(1:end-3));
matlabbatch{1}.spm.util.voi.roi{2}.mask.threshold = 0.0;
matlabbatch{1}.spm.util.voi.roi{2}.mask.image = cellstr(mask2);
matlabbatch{1}.spm.util.voi.expression        = 'i1 & i2';

% Run batch
spm_jobman('run',matlabbatch);

cd(start_dir);

% gzip
unix(['gzip ', fullfile(p1, p2)]);
