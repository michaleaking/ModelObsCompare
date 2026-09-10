function D = moc_spatial_diff(M, O, roi, varargin)
%MOC_SPATIAL_DIFF  Model-minus-observation speed difference on a common ROI grid.
%
%   Regrids the model mesh field and the gridded observations onto the SAME
%   regular grid over an ROI, then differences them. Handles the (common) case
%   where model and observations have different native grids and time sampling:
%   the observed frame(s) nearest a target time are averaged before comparison.
%
%   Usage:
%      D = moc_spatial_diff(M, O, roi)
%      D = moc_spatial_diff(M, O, roi, 'field','vel', 'res',200, ...
%                           'time',2019.0, 'obswindow',0.1)
%
%   Input:
%      M   : model struct (moc_load_model)
%      O   : obs grid struct (moc_load_obs_netcdf) with band 'vv' (or set 'obsband')
%      roi : ROI struct
%   Name/value options:
%      'field'     : (char)   model field ('vel' default). Compared to obs 'obsband'.
%      'obsband'   : (char)   obs band to compare (default 'vv' speed)
%      'res'       : (double) common-grid spacing [m] (default config_moc().grid_res)
%      'time'      : (double) target decimal year (default: model's last step)
%      'obswindow' : (double) +/- years of obs frames to average (default 0: nearest)
%
%   Output (struct D):
%      .x, .y     (1xnx / 1xny)   common grid
%      .X, .Y     (ny x nx)       meshgrid
%      .model     (ny x nx)       gridded model field [m/yr]
%      .obs       (ny x nx)       gridded observed field [m/yr]
%      .diff      (ny x nx)       model - obs [m/yr]
%      .time      (double)        model time used
%      .obstime   (1xk)           obs frame times averaged
%      .stats     (struct)        n, meanDiff, medianDiff, rmse, mae over valid pixels

p = inputParser;
p.addParameter('field','vel', @(s) ischar(s)||isstring(s));
p.addParameter('obsband','vv', @(s) ischar(s)||isstring(s));
p.addParameter('res',[], @(v) isempty(v)||isscalar(v));
p.addParameter('time',[], @(v) isempty(v)||isscalar(v));
p.addParameter('obswindow',0, @isscalar);
p.parse(varargin{:});
opt = p.Results;
obsband = char(opt.obsband);

% --- model on the common grid (also defines the grid) ---
G = moc_model_to_grid(M, roi, 'field', char(opt.field), 'res', opt.res, 'time', opt.time);

% --- select observed frame(s) near the target time ---
tTarget = opt.time;
if isempty(tTarget), tTarget = G.time; end
if isempty(tTarget) || isnan(tTarget), tTarget = median(O.time); end

if ~isfield(O, obsband)
    error('moc_spatial_diff:noBand', 'Obs struct has no band "%s".', obsband);
end
Ob = O.(obsband);   % (ny x nx x nt)

if opt.obswindow > 0
    sel = abs(O.time - tTarget) <= opt.obswindow;
    if ~any(sel), [~, k] = min(abs(O.time - tTarget)); sel = false(size(O.time)); sel(k)=true; end
else
    [~, k] = min(abs(O.time - tTarget));
    sel = false(size(O.time)); sel(k) = true;
end
obsSlice = mean(Ob(:,:,sel), 3, 'omitnan');    % (ny x nx) on the OBS grid
obstime = O.time(sel);

% --- interpolate observed slice onto the common (model) grid ---
% O.x is 1xnx (columns), O.y is 1xny (rows) -> obsSlice is (ny x nx).
[Ox, Oy] = meshgrid(O.x, O.y);
obsOnGrid = interp2(Ox, Oy, obsSlice, G.X, G.Y, 'linear', NaN);

% --- difference + stats over pixels valid in BOTH ---
diffField = G.grid - obsOnGrid;
valid = isfinite(diffField);
dv = diffField(valid);

D = struct();
D.x = G.x; D.y = G.y; D.X = G.X; D.Y = G.Y;
D.model = G.grid;
D.obs   = obsOnGrid;
D.diff  = diffField;
D.time  = G.time;
D.obstime = obstime;
D.stats = struct( ...
    'n',          numel(dv), ...
    'meanDiff',   mean(dv), ...
    'medianDiff', median(dv), ...
    'rmse',       sqrt(mean(dv.^2)), ...
    'mae',        mean(abs(dv)) );

fprintf(['Spatial diff (%s): model t=%.3f vs obs t=[%.3f..%.3f]\n' ...
         '  valid px=%d | mean=%.1f  median=%.1f  RMSE=%.1f  MAE=%.1f m/yr\n'], ...
    char(opt.field), D.time, min(obstime), max(obstime), ...
    D.stats.n, D.stats.meanDiff, D.stats.medianDiff, D.stats.rmse, D.stats.mae);

end
