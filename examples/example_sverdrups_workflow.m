%% example_sverdrups_workflow.m
% End-to-end tour of the toolkit on data that is actually on disk.
%
% Sverdrups Glacier (NW Greenland) is the one place in the configured folders
% where a gridded observation subset and a redrawn flowline cover the same ice:
%   obs grid : SentinelMonthly_Sverdrups.nc  (72 monthly frames, 2015-2020)
%   flowline : sverdrups.mat                 (31 points, 0-7.5 km, 2015-2024)
%
% Demonstrates, in order:
%   1. the two observation loaders,
%   2. building an ROI from the flowline footprint, saving it (.mat + ISSM .exp)
%      and reading it back,
%   3. a speed map with the ROI and flowline overlaid,
%   4. moc_velocity_variance: melt-season window, two individual time slices,
%      and a robust record-long variant + coverage map,
%   5. a per-year table of relative change over the ROI,
%   6. the flowline: distance-time image and a point time series.
% The model section (7) is OFF by default: it loads a multi-GB ISSM .mat.
%
% Usage:
%   >> example_sverdrups_workflow
%   >> SAVE_FIGS = true;  example_sverdrups_workflow    % also write figs/*.png
%   >> RUN_MODEL = true;  example_sverdrups_workflow    % add the model panels

%% 0. Options ---------------------------------------------------------------
% Pre-set any of these before running to override the defaults.
if ~exist('OBS_NC','var'),     OBS_NC     = 'SentinelMonthly_Sverdrups.nc'; end
if ~exist('OBS_FL','var'),     OBS_FL     = 'sverdrups.mat';                end
if ~exist('SAVE_FIGS','var'),  SAVE_FIGS  = false;                          end
if ~exist('RUN_MODEL','var'),  RUN_MODEL  = false;                          end
if ~exist('MODEL_FILE','var'), MODEL_FILE = 'Model_NW_HindcastRun_TransientInversion.mat'; end

cfg = setup_moc();
figdir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'figs');
if SAVE_FIGS && ~isfolder(figdir), mkdir(figdir); end

%% 1. Load the observations -------------------------------------------------
O = moc_load_obs_netcdf(OBS_NC, 'bands',{'vv'});
F = moc_load_obs_flowline(OBS_FL);

fprintf('\nobs grid  : %dx%d px, %d frames %.2f-%.2f, %.0f%% of pixel-frames valid\n', ...
    numel(O.x), numel(O.y), numel(O.time), min(O.time), max(O.time), ...
    100*mean(isfinite(O.vv(:))));
fprintf('flowline  : %s, %d points (0-%.1f km), %d times %.2f-%.2f (source %s)\n', ...
    F.name, numel(F.d), max(F.d)/1e3, numel(F.time), min(F.time), max(F.time), F.source);

%% 2. ROI from the flowline footprint, saved and read back ------------------
pad = 3000;                                   % [m] around the flowline
xr = [min(F.x)-pad, max(F.x)+pad];
yr = [min(F.y)-pad, max(F.y)+pad];
roi = struct('name','sverdrups_trunk', 'kind','rectangle', 'epsg',cfg.epsg, ...
             'x',[xr(1) xr(2) xr(2) xr(1) xr(1)], ...
             'y',[yr(1) yr(1) yr(2) yr(2) yr(1)]);
moc_save_roi(roi);                            % rois/sverdrups_trunk.mat + .exp
roi = moc_load_roi('sverdrups_trunk');        % round-trip through the saved file

[Xo, Yo] = meshgrid(O.x, O.y);
fprintf('ROI       : %.1f x %.1f km, %d of %d grid pixels inside\n', ...
    diff(xr)/1e3, diff(yr)/1e3, sum(moc_roi_mask(roi, Xo, Yo), 'all'), numel(Xo));

%% 3. A speed map with the ROI and flowline on top --------------------------
[~, ti] = min(abs(O.time - 2019.6));          % late-summer frame
f1 = figure('Color','w','Position',[80 80 720 600]);
ax = moc_plot_map(O, 'band','vv', 'tindex',ti, 'roi',roi, 'clim',[0 3000], ...
                  'ax',axes());   % without 'ax' the plotters open their own figure
plot(ax, F.x, F.y, 'w-', 'LineWidth',2);
plot(ax, F.x, F.y, 'k--', 'LineWidth',1);
if SAVE_FIGS, exportgraphics(f1, fullfile(figdir,'01_speed_map.png'), 'Resolution',120); end

%% 4a. Melt-season variability, normalised ----------------------------------
% Peak-to-peak speed over May-Sep 2019, divided by each pixel's own mean speed
% over the same window: the left panel is comparable between neighbouring
% pixels, the right one is dominated by the fast trunk simply because it is fast.
% 'minspeed' matters here — drop it and the ROI mean jumps from ~0.35 to ~1.2 as
% near-stagnant margin pixels divide a few m/yr of noise by a few m/yr of speed.
V = moc_velocity_variance(O, roi, 'trange',[2019.35 2019.85], ...
                          'mincount',3, 'minspeed',20);

f2 = figure('Color','w','Position',[80 80 1200 520]);
moc_plot_variance(V, 'field','norm',   'ax',subplot(1,2,1), 'roi',roi, ...
                  'title','relative: range / mean speed [-]');
moc_plot_variance(V, 'field','spread', 'ax',subplot(1,2,2), 'roi',roi, ...
                  'title','absolute: range of speed [m/yr]');
if SAVE_FIGS, exportgraphics(f2, fullfile(figdir,'02_melt_season_variability.png'), 'Resolution',120); end

%% 4b. Between two individual time slices -----------------------------------
% Nearest frame to each date; this window also yields a SIGNED change field.
V2 = moc_velocity_variance(O, roi, 'times',[2019.15 2019.6]);

f3 = figure('Color','w','Position',[80 80 1200 520]);
moc_plot_variance(V2, 'field','normChange', 'ax',subplot(1,2,1), 'roi',roi, ...
                  'title',sprintf('relative change %.2f -> %.2f [-]', V2.time));
moc_plot_variance(V2, 'field','change',     'ax',subplot(1,2,2), 'roi',roi, ...
                  'title',sprintf('absolute change %.2f -> %.2f [m/yr]', V2.time));
if SAVE_FIGS, exportgraphics(f3, fullfile(figdir,'03_two_slice_change.png'), 'Resolution',120); end
fprintf('\nlate winter (%.2f) -> late summer (%.2f): mean |change| = %.0f%% of local speed\n', ...
    V2.time(1), V2.time(2), 100*V2.stats.normMean);

%% 4c. Robust, record-long variant + coverage -------------------------------
% IQR about the median is insensitive to a single bad frame; 'minspeed' keeps
% near-stagnant ice from producing enormous ratios, 'mincount' demands enough
% frames per pixel. The count map shows where that bites.
Vr = moc_velocity_variance(O, roi, 'statistic','iqr', 'normalize','median', ...
                           'minspeed',20, 'mincount',24);

f4 = figure('Color','w','Position',[80 80 1200 520]);
moc_plot_variance(Vr, 'field','norm',  'ax',subplot(1,2,1), 'roi',roi, ...
                  'title',sprintf('2015-2020 IQR / median speed [-]  (mean %.2f)', ...
                                  Vr.stats.normMean));
moc_plot_variance(Vr, 'field','count', 'ax',subplot(1,2,2), 'roi',roi, ...
                  'title','valid frames per pixel (of 72)');
if SAVE_FIGS, exportgraphics(f4, fullfile(figdir,'04_record_long_robust.png'), 'Resolution',120); end

%% 5. Per-year relative change over the ROI ---------------------------------
% evalc() just swallows the per-call summary line so the table stays readable.
years = 2015:2020;
fprintf('\n  melt season (Apr-Nov) over the ROI\n');
fprintf('  year | frames | px    | mean rel. | median rel. | mean range [m/yr]\n');
fprintf('  -----+--------+-------+-----------+-------------+------------------\n');
for yr = years
    try
        evalc(sprintf(['Vy = moc_velocity_variance(O, roi, ''trange'',[%f %f], ' ...
                       '''mincount'',3, ''minspeed'',20);'], yr+0.25, yr+0.85));
        fprintf('  %4d |  %5d | %5d |   %6.1f%% |     %6.1f%% |   %8.0f\n', ...
            yr, numel(Vy.time), Vy.stats.npix, 100*Vy.stats.normMean, ...
            100*Vy.stats.normMedian, Vy.stats.spreadMean);
    catch ME
        fprintf('  %4d |  skipped (%s)\n', yr, ME.message);
    end
end

%% 6. The flowline ----------------------------------------------------------
% Distance-time image of observed speed, and the series at one point. (The
% model-vs-obs versions of these are moc_plot_flowline / moc_plot_timeseries,
% exercised in section 7.)
f5 = figure('Color','w','Position',[80 80 1250 480]);

ax1 = subplot(1,2,1);
imagesc(ax1, F.time, F.d/1e3, F.vel, 'AlphaData', ~isnan(F.vel));
set(ax1,'YDir','normal'); axis(ax1,'tight');
colormap(ax1,'parula'); cb = colorbar(ax1); ylabel(cb,'speed [m/yr]');
xlabel(ax1,'decimal year'); ylabel(ax1,'along-flowline distance [km]');
title(ax1, sprintf('%s flowline (%s)', F.name, F.source), 'Interpreter','none');

[~, kp] = min(abs(F.d - 3000));               % ~3 km up-glacier
ax2 = subplot(1,2,2); hold(ax2,'on'); grid(ax2,'on');
if ~isempty(F.err)
    errorbar(ax2, F.time, F.vel(kp,:), F.err(kp,:), 'o', 'MarkerSize',3, ...
             'Color',[0.2 0.2 0.2], 'CapSize',0);
else
    plot(ax2, F.time, F.vel(kp,:), 'o-', 'MarkerSize',3, 'Color',[0.2 0.2 0.2]);
end
% The same pixel from the gridded product, for an independent check.
[~, cx] = min(abs(O.x - F.x(kp)));
[~, cy] = min(abs(O.y - F.y(kp)));
plot(ax2, O.time, squeeze(O.vv(cy,cx,:)), 's-', 'MarkerSize',3, ...
     'Color',[0.85 0.1 0.1], 'LineWidth',1);
legend(ax2, {'flowline file','gridded obs (nearest px)'}, 'Location','best');
xlabel(ax2,'decimal year'); ylabel(ax2,'speed [m/yr]');
title(ax2, sprintf('%s at d=%.1f km', F.name, F.d(kp)/1e3), 'Interpreter','none');
if SAVE_FIGS, exportgraphics(f5, fullfile(figdir,'05_flowline.png'), 'Resolution',120); end

%% 7. Model comparison (optional — loads a multi-GB ISSM .mat) --------------
if RUN_MODEL
    M = moc_load_model(MODEL_FILE, 'solution','transient', 'trange',[2018 2020]);

    D = moc_spatial_diff(M, O, roi, 'field','vel', 'obsband','vv', ...
                         'res',cfg.grid_res, 'time',2019.6, 'obswindow',0.05);
    f6 = moc_plot_diff(D, 'roi',roi);
    if SAVE_FIGS, exportgraphics(f6, fullfile(figdir,'06_model_obs_diff.png'), 'Resolution',120); end

    Fmod = moc_model_flowline(M, F);
    f7 = figure('Color','w','Position',[80 80 700 520]);
    moc_plot_flowline(Fmod, F, 'mode','profile', 'time',2019.6, 'ax',axes());
    if SAVE_FIGS, exportgraphics(f7, fullfile(figdir,'07_flowline_profile.png'), 'Resolution',120); end

    C = moc_compare_timeseries(M, F, 'dist',3000);
    f8 = figure('Color','w','Position',[80 80 800 480]);
    moc_plot_timeseries(C, 'ax',axes());
    if SAVE_FIGS, exportgraphics(f8, fullfile(figdir,'08_point_timeseries.png'), 'Resolution',120); end
else
    fprintf(['\n(model section skipped — set RUN_MODEL = true to load %s\n' ...
             ' and add the model/obs difference, flowline profile and paired series)\n'], ...
            MODEL_FILE);
end

if SAVE_FIGS
    fprintf('\nfigures written to %s\n', figdir);
end
