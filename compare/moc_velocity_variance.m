function V = moc_velocity_variance(O, roi, varargin)
%MOC_VELOCITY_VARIANCE  Normalised temporal variability of speed, per pixel.
%
%   Answers "where, inside this ROI, did the ice speed change the most RELATIVE
%   to how fast it normally goes". The absolute spread (e.g. peak-to-peak speed
%   over a window) is dominated by the fast trunk of a glacier simply because it
%   is fast; dividing each pixel's spread by its own reference speed gives a
%   dimensionless map in which a pixel can be compared with its neighbours as a
%   fraction of relative change rather than an absolute magnitude change.
%
%   The temporal window is chosen in one of three ways:
%      'times'  : the frames nearest those decimal years (pass two values for
%                 the "two time slices" case; a signed change field is added)
%      'trange' : every frame inside that decimal-year window
%      neither  : every frame in O
%
%   Usage:
%      V = moc_velocity_variance(O)                        % full grid, all frames
%      V = moc_velocity_variance(O, roi, 'trange',[2019 2019.5])
%      V = moc_velocity_variance(O, roi, 'times',[2019.0 2019.5])
%      V = moc_velocity_variance(O, [], 'statistic','std', 'normalize','median')
%
%   Input:
%      O   : obs grid struct (moc_load_obs_netcdf) with band 'vv' (or set 'band')
%      roi : ROI struct, or [] for the full grid (may be omitted entirely)
%   Name/value options:
%      'band'      : (char)   band to use (default 'vv' speed)
%      'trange'    : (1x2)    [t0 t1] decimal-year window, inclusive
%      'times'     : (1xk)    target decimal years; nearest frame to each is used
%      'statistic' : (char)   spread measure over the window:
%                             'range' (default, max-min) | 'std' | 'iqr' | 'mad'
%      'normalize' : (char)   per-pixel reference speed:
%                             'mean' (default) | 'median' | 'first' | 'max' | 'none'
%      'mincount'  : (double) minimum finite frames a pixel needs (default 2)
%      'minspeed'  : (double) references below this [m/yr] are treated as invalid
%                             (default 0 = keep all). Stops near-stagnant ice
%                             from producing huge ratios.
%
%   Output (struct V):
%      .x, .y        (1xnx / 1xny)   ROI grid
%      .X, .Y        (ny x nx)       meshgrid
%      .spread       (ny x nx)       absolute spread [m/yr]
%      .reference    (ny x nx)       per-pixel normaliser [m/yr]
%      .norm         (ny x nx)       spread ./ reference [-] (the comparable map)
%      .count        (ny x nx)       finite frames per pixel
%      .change       (ny x nx)       later - earlier [m/yr]   } only when exactly
%      .normChange   (ny x nx)       change ./ reference [-]  } two frames are used
%      .time         (1xk)           frames used [decimal year]
%      .band, .statistic, .normalize (char), .roi (char), .epsg (double)
%      .stats        (struct)        npix, normMean/normMedian/normP90,
%                                    spreadMean/spreadMedian/spreadP90
%
%   See also MOC_PLOT_VARIANCE, MOC_SPATIAL_DIFF, MOC_ROI_MASK.

if nargin < 2, roi = []; end
if ischar(roi) || isstring(roi)      % roi omitted: moc_velocity_variance(O,'trange',...)
    varargin = [{roi}, varargin];
    roi = [];
end

STATISTICS = {'range','std','iqr','mad'};
NORMALIZERS = {'mean','median','first','max','none'};

p = inputParser;
p.addParameter('band','vv', @(s) ischar(s)||isstring(s));
p.addParameter('trange',[], @(v) isempty(v)||numel(v)==2);
p.addParameter('times',[], @(v) isempty(v)||isnumeric(v));
p.addParameter('statistic','range', @(s) any(strcmpi(char(s), STATISTICS)));
p.addParameter('normalize','mean', @(s) any(strcmpi(char(s), NORMALIZERS)));
p.addParameter('mincount',2, @isscalar);
p.addParameter('minspeed',0, @isscalar);
p.parse(varargin{:});
opt = p.Results;

band = char(opt.band);
statistic = lower(char(opt.statistic));
normalize = lower(char(opt.normalize));

if ~isfield(O, band)
    error('moc_velocity_variance:noBand', 'Obs struct has no band "%s".', band);
end

% --- ROI bounding-box crop (no-op for a full grid) -------------------------
x = O.x(:)';  y = O.y(:)';
if isempty(roi)
    xi = true(size(x));  yi = true(size(y));
    roiname = '';
else
    xi = x >= min(roi.x) & x <= max(roi.x);
    yi = y >= min(roi.y) & y <= max(roi.y);
    if ~any(xi) || ~any(yi)
        error('moc_velocity_variance:noOverlap', ...
              'ROI bbox (%.0f..%.0f, %.0f..%.0f) does not overlap the grid.', ...
              min(roi.x), max(roi.x), min(roi.y), max(roi.y));
    end
    if isfield(roi,'name'), roiname = char(roi.name); else, roiname = ''; end
end
x = x(xi);  y = y(yi);

% --- temporal window -> frame indices (time-ascending) --------------------
tAll = O.time(:)';
k = local_time_indices(tAll, opt.trange, opt.times);
A = double(O.(band)(yi, xi, k));       % (ny x nx x nt) over the ROI window
t = tAll(k);

% --- per-pixel spread and reference ---------------------------------------
count = sum(isfinite(A), 3);
keep  = count >= max(round(opt.mincount), 1);

spread    = local_spread(A, statistic);
reference = local_reference(A, normalize);

badRef = ~isfinite(reference) | reference <= 0 | reference < opt.minspeed;
if strcmp(normalize, 'none')
    nrm = spread;
else
    nrm = spread ./ reference;
    nrm(badRef) = NaN;
end

% --- polygon mask (a rectangle ROI is already handled by the bbox crop) ----
[X, Y] = meshgrid(x, y);
if ~isempty(roi)
    keep = keep & moc_roi_mask(roi, X, Y);
end

spread(~keep)    = NaN;
reference(~keep) = NaN;
nrm(~keep)       = NaN;

V = struct();
V.x = x;  V.y = y;  V.X = X;  V.Y = Y;
V.spread    = spread;
V.reference = reference;
V.norm      = nrm;
V.count     = count;

if numel(t) == 2
    change = A(:,:,2) - A(:,:,1);       % later minus earlier
    change(~keep) = NaN;
    normChange = change ./ reference;
    normChange(badRef | ~keep) = NaN;
    V.change     = change;
    V.normChange = normChange;
end

V.time      = t;
V.band      = band;
V.statistic = statistic;
V.normalize = normalize;
V.mincount  = opt.mincount;
V.minspeed  = opt.minspeed;
V.roi       = roiname;
if isfield(O,'epsg'), V.epsg = O.epsg; else, V.epsg = NaN; end

sN = local_summary(nrm);
sS = local_summary(spread);
V.stats = struct( ...
    'npix',         sN.npix, ...
    'normMean',     sN.mean,   'normMedian',   sN.median,   'normP90',   sN.p90, ...
    'spreadMean',   sS.mean,   'spreadMedian', sS.median,   'spreadP90', sS.p90);

fprintf(['Velocity variance (%s of %s, normalised by %s): %d frame(s), ' ...
         't=[%.3f..%.3f]\n' ...
         '  valid px=%d | norm mean=%.3f  median=%.3f  p90=%.3f  ' ...
         '| spread mean=%.1f m/yr\n'], ...
    statistic, band, normalize, numel(t), min(t), max(t), ...
    V.stats.npix, V.stats.normMean, V.stats.normMedian, V.stats.normP90, ...
    V.stats.spreadMean);

end

% ===========================================================================
function k = local_time_indices(t, trange, times)
%LOCAL_TIME_INDICES  Resolve the window options to sorted, unique frame indices.
if ~isempty(times) && ~isempty(trange)
    error('moc_velocity_variance:bothWindows', ...
          'Pass either ''times'' or ''trange'', not both.');
end

if ~isempty(times)
    targets = times(:)';
    k = zeros(1, numel(targets));
    for i = 1:numel(targets)
        [~, k(i)] = min(abs(t - targets(i)));
    end
    k = unique(k);
    if numel(k) < numel(targets)
        warning('moc_velocity_variance:collapsedTimes', ...
            ['%d target times collapsed onto %d distinct frame(s); the ' ...
             'requested slices are not resolved by this series.'], ...
            numel(targets), numel(k));
    end
elseif ~isempty(trange)
    k = find(t >= min(trange) & t <= max(trange));
    if isempty(k)
        error('moc_velocity_variance:noFrames', ...
              'No frames in trange [%g %g]; series spans %.3f..%.3f.', ...
              min(trange), max(trange), min(t), max(t));
    end
else
    k = 1:numel(t);
end

if numel(k) < 2
    error('moc_velocity_variance:tooFewFrames', ...
          'Need at least 2 time slices to measure variability.');
end
[~, ord] = sort(t(k));   % time-ascending, so 'first'/change are ordered
k = k(ord);
end

% ---------------------------------------------------------------------------
function s = local_spread(A, statistic)
%LOCAL_SPREAD  Reduce a (ny x nx x nt) stack to a per-pixel spread [m/yr].
switch statistic
    case 'range'
        s = max(A, [], 3, 'omitnan') - min(A, [], 3, 'omitnan');
    case 'std'
        s = std(A, 0, 3, 'omitnan');
    case 'iqr'
        s = local_pctl3(A, 75) - local_pctl3(A, 25);
    case 'mad'
        s = median(abs(A - median(A, 3, 'omitnan')), 3, 'omitnan');
    otherwise
        error('moc_velocity_variance:badStatistic', 'Unknown statistic "%s".', statistic);
end
s(all(~isfinite(A), 3)) = NaN;   % all-NaN pixels: max-min would give -Inf
end

% ---------------------------------------------------------------------------
function r = local_reference(A, normalize)
%LOCAL_REFERENCE  Per-pixel reference speed [m/yr] used to normalise the spread.
%   'none' still returns the mean, reported for context but not divided by.
switch normalize
    case 'median'
        r = median(A, 3, 'omitnan');
    case 'max'
        r = max(A, [], 3, 'omitnan');
    case 'first'
        % Earliest FINITE frame per pixel (individual frames may be NaN).
        [ny, nx, ~] = size(A);
        [~, first] = max(isfinite(A), [], 3);     % max returns the first true
        r = NaN(ny, nx);
        idx = (1:ny*nx)';
        r(idx) = A(idx + (first(idx)-1)*ny*nx);
    otherwise   % 'mean' and 'none'
        r = mean(A, 3, 'omitnan');
end
r(all(~isfinite(A), 3)) = NaN;
end

% ---------------------------------------------------------------------------
function q = local_pctl3(A, pc)
%LOCAL_PCTL3  Nearest-rank percentile along dim 3, ignoring NaNs (no toolbox).
[ny, nx, ~] = size(A);
S = sort(A, 3);                      % NaNs sort to the end
n = sum(isfinite(A), 3);
q = NaN(ny, nx);
good = find(n > 0);
rank = max(1, ceil(pc/100 * n(good)));
q(good) = S(good + (rank-1)*ny*nx);
end

% ---------------------------------------------------------------------------
function s = local_summary(F)
%LOCAL_SUMMARY  Mean / median / 90th-percentile of a 2-D field's finite pixels.
v = F(isfinite(F));
if isempty(v)
    s = struct('npix',0,'mean',NaN,'median',NaN,'p90',NaN);
    return
end
sv = sort(v);
s = struct('npix', numel(v), 'mean', mean(v), 'median', median(v), ...
           'p90', sv(max(1, round(0.90*numel(sv)))));
end
