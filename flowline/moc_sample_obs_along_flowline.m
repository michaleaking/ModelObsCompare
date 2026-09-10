function Fobs = moc_sample_obs_along_flowline(O, FL, varargin)
%MOC_SAMPLE_OBS_ALONG_FLOWLINE  Sample a gridded observation struct along a flowline.
%
%   Turns the geometry-only struct from MOC_FLOWLINE_FROM_ROI into the same
%   layout MOC_LOAD_OBS_FLOWLINE returns, so it can be handed to
%   MOC_PLOT_FLOWLINE and MOC_COMPARE_TIMESERIES.
%
%   Usage:
%      Fobs = moc_sample_obs_along_flowline(O, FL)
%      Fobs = moc_sample_obs_along_flowline(O, FL, 'band','vv', 'method','nearest')
%
%   Input:
%      O  : gridded obs struct from MOC_LOAD_OBS_NETCDF (.x, .y, .time, .(band))
%      FL : flowline struct with .x, .y, .d
%   Name/value options:
%      'band'   : (char) band to sample as .vel (default 'vv'); vx/vy/ex/ey are
%                 carried across too when present in O
%      'method' : (char) interp2 method, 'linear'(default) | 'nearest'
%
%   Output (struct Fobs): as MOC_LOAD_OBS_FLOWLINE —
%      .x, .y, .d  (np x 1)     the flowline
%      .time       (1 x nt)     decimal years
%      .vel        (np x nt)    sampled speed [m/yr]
%      .vx,.vy,.err(np x nt)    where available ([] otherwise)
%      .source     (char)       the band sampled
%      .name, .file, .epsg

p = inputParser;
p.addParameter('band', 'vv', @(s) ischar(s)||isstring(s));
p.addParameter('method', 'linear', @(s) ischar(s)||isstring(s));
p.parse(varargin{:});
band = char(p.Results.band);
method = char(p.Results.method);

if ~isfield(O, band)
    error('moc_sample_obs_along_flowline:noBand', ...
        'Band "%s" not in the obs struct (have: %s).', band, ...
        strjoin(intersect(fieldnames(O)', {'vv','vx','vy','ex','ey','dT'}), ', '));
end

px = FL.x(:); py = FL.y(:);
[Ox, Oy] = meshgrid(O.x, O.y);

Fobs = struct();
Fobs.name = ''; if isfield(FL,'name'), Fobs.name = FL.name; end
Fobs.file = ''; if isfield(O,'file'),  Fobs.file = O.file;  end
Fobs.x = px; Fobs.y = py; Fobs.d = FL.d(:);
Fobs.time = O.time(:)';
Fobs.source = band;
Fobs.epsg = FL.epsg;

Fobs.vel = local_sample(O.(band), Ox, Oy, px, py, method);
for b = {'vx','vy'}
    if isfield(O, b{1}) && ~strcmp(b{1}, band)
        Fobs.(b{1}) = local_sample(O.(b{1}), Ox, Oy, px, py, method);
    else
        Fobs.(b{1}) = [];
    end
end
Fobs.err = [];
if isfield(O, 'ex') && isfield(O, 'ey')
    ex = local_sample(O.ex, Ox, Oy, px, py, method);
    ey = local_sample(O.ey, Ox, Oy, px, py, method);
    Fobs.err = hypot(ex, ey);
end

% Carry the provenance of the clip across
for f = {'glacier','flowline','terminus_date','terminus_file'}
    if isfield(FL, f{1}), Fobs.(f{1}) = FL.(f{1}); end
end

end

% ===========================================================================
function V = local_sample(B, Ox, Oy, px, py, method)
%LOCAL_SAMPLE  Interpolate an (ny x nx x nt) band onto np points -> (np x nt).
nt = size(B, 3);
V = nan(numel(px), nt);
for k = 1:nt
    V(:,k) = interp2(Ox, Oy, B(:,:,k), px, py, method, NaN);
end
end
