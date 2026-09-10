function C = moc_compare_timeseries(M, obs, varargin)
%MOC_COMPARE_TIMESERIES  Align modelled and observed velocity time series at a point.
%
%   Extracts a modelled time series at a chosen location and pairs it with an
%   observed series taken either from a flowline struct (nearest along-line
%   point) or from a gridded obs struct (nearest pixel). Returns both series
%   plus the model interpolated onto the observation times for error metrics.
%
%   Usage:
%      % against an observed flowline, at a given along-flowline distance:
%      C = moc_compare_timeseries(M, Fobs, 'dist',5000)
%      % against an observed flowline, at nearest point to a coordinate:
%      C = moc_compare_timeseries(M, Fobs, 'xy',[x y])
%      % against a gridded obs struct, nearest pixel to a coordinate:
%      C = moc_compare_timeseries(M, O, 'xy',[x y], 'obsband','vv')
%
%   Input:
%      M   : model struct (moc_load_model)
%      obs : observed flowline struct (moc_load_obs_flowline) OR
%            gridded obs struct (moc_load_obs_netcdf)
%   Name/value options:
%      'xy'      : (1x2) [x y] location [m] (required for grid obs; optional for flowline)
%      'dist'    : (double) along-flowline distance [m] to sample (flowline obs only)
%      'field'   : (char) model field ('vel' default)
%      'obsband' : (char) obs band for grid obs ('vv' default)
%
%   Output (struct C):
%      .xy            (1x2)   location used [m]
%      .model_time    (1xnt)  decimal years (model)
%      .model_val     (1xnt)  modelled series
%      .obs_time      (1xno)  decimal years (obs)
%      .obs_val       (1xno)  observed series
%      .obs_err       (1xno)  observed error ([] if unavailable)
%      .model_at_obs  (1xno)  model interpolated onto obs times
%      .stats         struct  bias, rmse, mae, r (over co-valid obs times)

p = inputParser;
p.addParameter('xy',[], @(v) isempty(v)||numel(v)==2);
p.addParameter('dist',[], @(v) isempty(v)||isscalar(v));
p.addParameter('field','vel', @(s) ischar(s)||isstring(s));
p.addParameter('obsband','vv', @(s) ischar(s)||isstring(s));
p.parse(varargin{:});
opt = p.Results;
field = char(opt.field);

isFlow = isstruct(obs) && isfield(obs,'d') && isfield(obs,'vel') && isfield(obs,'source');

% ---- pick location + observed series ----
if isFlow
    [ix, xy, obs_time, obs_val, obs_err] = local_flowline_pick(obs, opt);
else
    [xy, obs_time, obs_val, obs_err] = local_grid_pick(obs, opt, char(opt.obsband));
end

% ---- modelled series at that location ----
Pm = moc_model_at_points(M, xy(1), xy(2), 'field', field);
model_time = Pm.time;
model_val  = Pm.val(1, :);

% ---- model on obs times + stats ----
if numel(model_time) >= 2 && ~any(isnan(model_time))
    model_at_obs = interp1(model_time, model_val, obs_time, 'linear', NaN);
else
    model_at_obs = nan(size(obs_time));
end
good = isfinite(model_at_obs) & isfinite(obs_val);
res = model_at_obs(good) - obs_val(good);

C = struct();
C.xy = xy;
C.model_time = model_time;
C.model_val = model_val;
C.obs_time = obs_time;
C.obs_val = obs_val;
C.obs_err = obs_err;
C.model_at_obs = model_at_obs;
if any(good)
    r = corrcoef(model_at_obs(good), obs_val(good));
    C.stats = struct('n',nnz(good), 'bias',mean(res), 'rmse',sqrt(mean(res.^2)), ...
                     'mae',mean(abs(res)), 'r', r(1,2));
else
    C.stats = struct('n',0,'bias',NaN,'rmse',NaN,'mae',NaN,'r',NaN);
end

end

% ===========================================================================
function [ix, xy, t, v, e] = local_flowline_pick(F, opt)
%LOCAL_FLOWLINE_PICK  Choose a flowline point by distance or nearest coordinate.
if ~isempty(opt.dist)
    [~, ix] = min(abs(F.d - opt.dist));
elseif ~isempty(opt.xy)
    [~, ix] = min(hypot(F.x - opt.xy(1), F.y - opt.xy(2)));
else
    ix = 1;   % default: terminus / first point
end
xy = [F.x(ix), F.y(ix)];
t  = F.time(:)';
v  = F.vel(ix, :);
e  = [];
if isfield(F,'err') && ~isempty(F.err), e = F.err(ix, :); end
end

% ---------------------------------------------------------------------------
function [xy, t, v, e] = local_grid_pick(O, opt, band)
%LOCAL_GRID_PICK  Nearest-pixel time series from a gridded obs struct.
if isempty(opt.xy)
    error('moc_compare_timeseries:needXY', ...
        'For gridded obs you must give ''xy'',[x y].');
end
[~, ci] = min(abs(O.x - opt.xy(1)));
[~, ri] = min(abs(O.y - opt.xy(2)));
xy = [O.x(ci), O.y(ri)];
if ~isfield(O, band)
    error('moc_compare_timeseries:noBand', 'Obs struct has no band "%s".', band);
end
v = squeeze(O.(band)(ri, ci, :))';
t = O.time(:)';
e = [];
if isfield(O,'ex') && isfield(O,'ey')
    e = squeeze(hypot(O.ex(ri,ci,:), O.ey(ri,ci,:)))';
end
end
