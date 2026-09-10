%EXAMPLE_FLOWLINE_FROM_SHAPEFILES  Build a comparison-ready flowline from the shapefiles.
%
%   Takes a region of interest and produces the single centre-glacier flowline
%   that runs through it, clipped at the most recently traced terminus position:
%
%      ROI  ->  Felikson flowlines that reach into it
%           ->  centre flowline of the busiest glacier
%           ->  split at the newest Black & Joughin terminus trace, upstream side
%           ->  clipped to the ROI, resampled
%           ->  flowline struct -> moc_model_flowline / moc_plot_flowline
%
%   MATLAB mirror of examples/example_flowline_from_shapefiles.py.
%   Needs the Mapping Toolbox (shaperead, polyxpoly).
%
%   Usage:
%      >> setup_moc
%      >> example_flowline_from_shapefiles
%      >> ROI_NAME = 'my_roi'; example_flowline_from_shapefiles
%      >> GLACIER = 'b045'; FLOWLINE = '05'; example_flowline_from_shapefiles
%      >> RUN_MODEL = true; example_flowline_from_shapefiles   % add the model panel
%
%   The model section is off by default: it loads a multi-GB ISSM .mat, where
%   every other section runs in a second or two.

if ~exist('ROI_NAME','var'),  ROI_NAME  = 'sverdrups_trunk';                 end
if ~exist('GLACIER','var'),   GLACIER   = '';                                end
if ~exist('FLOWLINE','var'),  FLOWLINE  = '';                                end
if ~exist('SPACING','var'),   SPACING   = 200;                               end
if ~exist('OBS_NC','var'),    OBS_NC    = 'SentinelMonthly_Sverdrups.nc';    end
if ~exist('RUN_MODEL','var'), RUN_MODEL = false;                             end
if ~exist('MODEL_FILE','var'),MODEL_FILE= 'Model_NW_HindcastRun_TransientInversion.mat'; end
PROFILE_TIME = 2019.6;

cfg = config_moc();
roi = moc_load_roi(ROI_NAME);

%% 1. what the ROI catches, before any clipping
F = moc_load_flowlines(roi);
T = moc_load_termini(roi);
fprintf('\nflowlines : %d in the ROI, glaciers %s\n', numel(F), ...
    strjoin(unique({F.glacier}), ', '));
for g = unique({F.glacier})
    ids = {F(strcmp({F.glacier}, g{1})).flowline};
    fprintf('    %s: flowlines %s\n', g{1}, strjoin(sort(ids), ', '));
end
if isempty(T)
    fprintf('termini   : none in the ROI\n');
else
    fprintf('termini   : %d traces, %s to %s\n', numel(T), ...
        datestr(min([T.date]),'yyyy-mm-dd'), datestr(max([T.date]),'yyyy-mm-dd')); %#ok<DATST>
end

%% 2. the pipeline
FL = moc_flowline_from_roi(roi, 'glacier', GLACIER, 'flowline', FLOWLINE, ...
                           'spacing', SPACING);
fprintf('\nflowline  : glacier %s line %s, %d points over %.2f km\n', ...
    FL.glacier, FL.flowline, numel(FL.x), FL.length/1e3);
fprintf('            clipped at the %s terminus trace\n', FL.terminus_date);
fprintf('            d = 0 at (%.0f, %.0f)\n', FL.x(1), FL.y(1));

%% 3. map the ROI, its flowlines, and the terminus used for the clip
[rx, ry] = moc_roi_polygon(roi);
figure('Color','w','Name','flowlines and termini in the ROI');
ax = axes(); hold(ax,'on');
hRoi = plot(ax, rx, ry, 'k-', 'LineWidth', 1.5);
for k = 1:numel(F)
    plot(ax, F(k).x, F(k).y, '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
end
for k = 1:numel(T)
    plot(ax, T(k).x, T(k).y, '-', 'Color', [0.50 0.70 0.90], 'LineWidth', 0.5);
end
hTerm = [];
if ~isempty(T)
    [tNew, ~] = moc_latest_terminus(T, 'line', [FL.x FL.y]);
    if isempty(tNew), [~, kNew] = max([T.date]); tNew = T(kNew); end
    hTerm = plot(ax, tNew.x, tNew.y, '-', 'Color', [0.12 0.44 0.71], 'LineWidth', 2.5);
end
hFl  = plot(ax, FL.x, FL.y, '-', 'Color', [0.85 0.11 0.11], 'LineWidth', 2.5);
hEnd = plot(ax, FL.x(1), FL.y(1), 'o', 'MarkerFaceColor', [0.85 0.11 0.11], ...
            'MarkerEdgeColor', 'none', 'MarkerSize', 7);
axis(ax, 'equal');
pad = 0.15 * max(range(rx), range(ry));
xlim(ax, [min(rx)-pad max(rx)+pad]);
ylim(ax, [min(ry)-pad max(ry)+pad]);
xlabel(ax, 'x [m, EPSG:3413]'); ylabel(ax, 'y [m, EPSG:3413]');
title(ax, 'flowlines and termini in the ROI');
if isempty(hTerm)
    legend(ax, [hRoi hFl hEnd], {'ROI', ...
        sprintf('clipped flowline (%s)', FL.name), 'd = 0 (terminus)'}, ...
        'Location','best','FontSize',8,'Interpreter','none');
else
    legend(ax, [hRoi hTerm hFl hEnd], {'ROI', ...
        sprintf('terminus %s', FL.terminus_date), ...
        sprintf('clipped flowline (%s)', FL.name), 'd = 0 (terminus)'}, ...
        'Location','best','FontSize',8,'Interpreter','none');
end
hold(ax,'off');

%% 4. hand it the observations
Fobs = [];
try
    O = moc_load_obs_netcdf(OBS_NC, 'bands', {'vv'});
catch ME
    fprintf('\n(no gridded obs at %s — stopping after the geometry: %s)\n', ...
        OBS_NC, ME.message);
    O = [];
end
if ~isempty(O)
    Fobs = moc_sample_obs_along_flowline(O, FL);
    fprintf('\nobs       : %d frames sampled onto the line, %.0f%% of point-frames valid\n', ...
        numel(Fobs.time), 100*mean(isfinite(Fobs.vel(:))));

    % observed distance-time image (moc_plot_flowline needs a model, so draw it here)
    figure('Color','w','Name','observed speed along the flowline');
    imagesc(Fobs.time, Fobs.d/1e3, Fobs.vel);
    set(gca,'YDir','normal'); colormap(gca, parula); cb = colorbar();
    cb.Label.String = 'speed [m/yr]';
    xlabel('time [decimal year]'); ylabel('along-flowline distance [km]');
    title(sprintf('%s: observed speed (%s)', FL.name, Fobs.source), 'Interpreter','none');
end

%% 5. and the model
if RUN_MODEL && ~isempty(Fobs)
    M = moc_load_model(MODEL_FILE, 'solution','transient', 'trange',[2018 2020]);
    fprintf('model     : %d vertices, %d steps\n', numel(M.x), numel(M.time));
    Fmod = moc_model_flowline(M, FL);
    C = moc_compare_timeseries(M, Fobs, 'dist', 0.5*FL.length);
    fprintf('model-obs : at d=%.1f km  bias=%.0f  rmse=%.0f m/yr  r=%.2f over %d paired times\n', ...
        0.5*FL.length/1e3, C.stats.bias, C.stats.rmse, C.stats.r, C.stats.n);
    figure('Color','w','Name','model vs obs along the flowline');
    moc_plot_flowline(Fmod, Fobs, 'mode','profile', 'time', PROFILE_TIME);
elseif ~RUN_MODEL
    fprintf('\n(model section skipped — set RUN_MODEL = true to load %s)\n', MODEL_FILE);
end
