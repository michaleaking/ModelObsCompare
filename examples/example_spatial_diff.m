%% example_spatial_diff.m
% Compare modelled vs observed surface speed as a 2-D difference map over an ROI.
% Pipeline: load model -> load gridded obs -> pick/define ROI -> regrid both to a
% common grid -> difference -> plot three panels (model / obs / model-obs).

cfg = setup_moc();

%% 1. Pick a model run ------------------------------------------------------
T = moc_list_models();          % table of available runs with parsed metadata
disp(T);

% Load a transient run (keep only a couple of years to save memory if desired).
M = moc_load_model('Model_NW_HindcastRun_SnapshotInversion2007_fric1.mat', ...
                   'solution','transient');   % add 'trange',[2018 2020] to subset

%% 2. Load gridded observations for the same area ---------------------------
O = moc_load_obs_netcdf('Sentinel_Subset_Upernavik.nc', ...
                        'bands',{'vv'}, 'trange',[2018 2020]);

%% 3. Define the ROI --------------------------------------------------------
% Reuse a saved ROI, or derive one from the obs extent:
if isfile(fullfile(moc_path('roi_dir'),'upernavik_box.mat'))
    roi = moc_load_roi('upernavik_box');
else
    roi = moc_roi_from_obs(O, 'margin',2000);
end

%% 4. Regrid + difference at a target time ----------------------------------
targetYear = 2019.6;
D = moc_spatial_diff(M, O, roi, 'field','vel', 'obsband','vv', ...
                     'res', cfg.grid_res, 'time', targetYear, 'obswindow', 0.1);

%% 5. Plot ------------------------------------------------------------------
moc_plot_diff(D, 'roi', roi);
% D.stats holds n / meanDiff / medianDiff / rmse / mae over co-valid pixels.
fprintf('RMSE=%.1f  MAE=%.1f  bias=%.1f m/yr over %d px\n', ...
    D.stats.rmse, D.stats.mae, D.stats.meanDiff, D.stats.n);
