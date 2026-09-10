%% example_flowline_timeseries.m
% Compare modelled vs observed velocity time series along a glacier flowline.
% Uses the per-glacier observed flowline .mat files (redrawn_Jan25), which give
% both the (x,y) flowline geometry AND the observed velocity time series.

cfg = setup_moc();

%% 1. Load an observed flowline --------------------------------------------
% Auto-picks the first available velocity source (filteredV / SentinelV / u2 ...).
Fobs = moc_load_obs_flowline('tracy.mat');   % or 'source','SentinelV'

%% 2. Load a model run ------------------------------------------------------
M = moc_load_model('Model_NW_HindcastRun_SnapshotInversion2007_fric1.mat', ...
                   'solution','transient');

%% 3. Sample the model along the SAME flowline ------------------------------
Fmod = moc_model_flowline(M, Fobs, 'field','vel');   % reuses Fobs.x, Fobs.y, Fobs.d

%% 4a. Profile comparison at one time --------------------------------------
moc_plot_flowline(Fmod, Fobs, 'mode','profile', 'time',2019.0);

%% 4b. Distance-time difference (Hovmoller) --------------------------------
moc_plot_flowline(Fmod, Fobs, 'mode','hovmoller');

%% 5. Single-point time series at a chosen along-flowline distance ----------
% e.g. 5 km up-glacier from the terminus
C = moc_compare_timeseries(M, Fobs, 'dist',3000, 'field','vel');
moc_plot_timeseries(C, 'datetime',true);
fprintf('At d=5 km: RMSE=%.0f  bias=%.0f  r=%.2f  (n=%d)\n', ...
    C.stats.rmse, C.stats.bias, C.stats.r, C.stats.n);

%% 6. Or compare against a gridded obs pixel instead of the flowline file ---
% O = moc_load_obs_netcdf('Sentinel_Subset_Upernavik.nc','bands',{'vv','ex','ey'});
% C2 = moc_compare_timeseries(M, O, 'xy',[Fobs.x(20) Fobs.y(20)], 'obsband','vv');
% moc_plot_timeseries(C2, 'datetime',true);
