%% example_export_model.m
% Phase 2 bridge: export an ISSM model's mesh + velocity to a neutral netCDF
% that the Python toolkit (moc_py) can read with h5py. Run this in MATLAB once
% per model; then do the comparison/plotting in a Jupyter notebook.

cfg = setup_moc();

%% 1. Load a model (subset transient time to keep the export small if desired)
% A bare filename is resolved against moc_path('model_dir'), so this runs
% wherever your models live. Override MODEL_FILE before running for another one,
% or use moc_list_models() to see what is in the configured folder.
if ~exist('MODEL_FILE','var')
    MODEL_FILE = 'Model_NW_HindcastRun_SnapshotInversion2007_fric1.mat';
end
M = moc_load_model(MODEL_FILE, 'solution','transient');%, 'trange',[2018 2020]);

%% 2. Export to netCDF (mesh x/y/elements, time, vel/vx/vy [+surface/bed])
outdir = fullfile(moc_path('model_dir'), 'exported');
if ~isfolder(outdir), mkdir(outdir); end
ncpath = fullfile(outdir, 'Model_NW_fric1_Snapshot.nc');

moc_export_model(M, ncpath, 'fields',{'vel','vx','vy'});

% Snapshot inversions work too (single time step, time=NaN):
% Ms = moc_load_model('Model_NW_fric1_SnapshotInversion_2007.mat','solution','stressbalance');
% moc_export_model(Ms, fullfile(outdir,'Model_NW_snap2007.nc'));

fprintf('\nNow in Python:\n');
fprintf('  import moc_py as mp\n');
fprintf('  M = mp.load_model(%s)\n', ['"' ncpath '"']);
fprintf('  D = mp.spatial_diff(M, O, roi, time=2019.6)\n');
