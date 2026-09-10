function F = moc_load_obs_flowline(matpath, varargin)
%MOC_LOAD_OBS_FLOWLINE  Load a per-glacier observed flowline .mat (redrawn_Jan25).
%
%   These files describe ONE glacier flowline: a set of (x,y) points ordered
%   from the terminus up-glacier, with an along-flowline distance `d`, static
%   geometry (bed `zb`, thickness `zthick`), and one or more observed velocity
%   time series stored as points-by-time (or time-by-points) matrices inside
%   structs. This loader normalises the common variants to a single layout.
%
%   Recognised velocity sources (first present wins unless 'source' given):
%      filteredV        : (npts x nt)          paired with `vti`  (datenum)
%      SentinelV        : struct .velocity (npts x nt), .err, .t (datenum)
%      SentinelVmonthly : struct .velocity/.vx/.vy (npts x nt), .t (datenum)
%      u2               : struct .v (nt x npts!), .ve, .t (datenum)  [auto-transposed]
%
%   Usage:
%      F = moc_load_obs_flowline('jakobshavn.mat')      % bare name -> obs_flowline_dir
%      F = moc_load_obs_flowline(matpath, 'source','SentinelV')
%
%   Name/value options:
%      'source' : (char) which velocity source to extract
%                 'auto'(default)|'filteredV'|'SentinelV'|'SentinelVmonthly'|'u2'
%
%   Output (struct F):
%      .file    (char)             source file
%      .name    (char)             glacier name
%      .x, .y   (npts x 1 double)  flowline coordinates [m], EPSG per config
%      .d       (npts x 1 double)  along-flowline distance [m] (0 at first point)
%      .bed     (npts x 1 double)  bed elevation [m]           (zb, if present)
%      .thick   (npts x 1 double)  ice thickness [m]           (zthick, if present)
%      .source  (char)             velocity source used
%      .time    (1 x nt double)    decimal years
%      .vel     (npts x nt double) observed speed [m/yr]
%      .err     (npts x nt double) speed error [m/yr] ([] if unavailable)
%      .vx,.vy  (npts x nt double) components [m/yr] ([] if unavailable)

p = inputParser;
p.addParameter('source', 'auto', @(s) ischar(s)||isstring(s));
p.parse(varargin{:});
srcReq = char(p.Results.source);

% Resolve path
if ~isfile(matpath)
    cfg = config_moc();
    cand = fullfile(moc_path('obs_flowline_dir'), matpath);
    if isfile(cand), matpath = cand; else
        error('moc_load_obs_flowline:notFound', 'Flowline file not found: %s', matpath);
    end
end

S = load(matpath);

F = struct();
F.file = matpath;
[~, base] = fileparts(matpath);
F.name = base;
if isfield(S,'glrname'), try, F.name = strtrim(char(S.glrname)); catch, end; end

F.x    = local_col(S, 'x');
F.y    = local_col(S, 'y');
F.d    = local_col(S, 'd');
F.bed  = local_col(S, 'zb');
F.thick= local_col(S, 'zthick');
npts   = numel(F.x);

% Choose a velocity source
srcOrder = {'filteredV','SentinelV','SentinelVmonthly','u2'};
if ~strcmpi(srcReq,'auto')
    srcOrder = {srcReq};
end

F.source = '';
F.vel = []; F.err = []; F.vx = []; F.vy = []; F.time = [];

for k = 1:numel(srcOrder)
    src = srcOrder{k};
    if ~isfield(S, src), continue; end
    [vel, err, vx, vy, tdn, ok] = local_extract(S, src, npts);
    if ok
        F.source = src;
        F.vel  = vel;
        F.err  = err;
        F.vx   = vx;
        F.vy   = vy;
        F.time = moc_datenum2decyear(tdn(:)');
        break
    end
end

if isempty(F.source)
    warning('moc_load_obs_flowline:noVel', ...
        'No usable velocity source found in %s (looked for %s).', matpath, strjoin(srcOrder,', '));
end

fprintf('Loaded flowline: %s | %d pts | source=%s | %d frames\n', ...
    F.name, npts, F.source, numel(F.time));

end

% ===========================================================================
function v = local_col(S, name)
%LOCAL_COL  Fetch S.(name) as a column vector, or [] if absent.
v = [];
if isfield(S, name)
    v = S.(name);
    v = v(:);
end
end

% ---------------------------------------------------------------------------
function [vel, err, vx, vy, tdn, ok] = local_extract(S, src, npts)
%LOCAL_EXTRACT  Pull (npts x nt) velocity/err/components + datenum time from a source.
vel = []; err = []; vx = []; vy = []; tdn = []; ok = false;

switch src
    case 'filteredV'
        vel = S.filteredV;                 % (npts x nt)
        if isfield(S,'filteredVe'), err = S.filteredVe; end
        if isfield(S,'vti'), tdn = S.vti(:); end

    case {'SentinelV','SentinelVmonthly'}
        st = S.(src);
        if isstruct(st) && numel(st)>=1, st = st(1); end
        vel = local_getf(st, 'velocity');
        err = local_getf(st, 'err');
        vx  = local_getf(st, 'vx');
        vy  = local_getf(st, 'vy');
        tdn = local_getf(st, 't');

    case 'u2'
        st = S.u2;
        if isstruct(st) && numel(st)>=1, st = st(1); end
        vel = local_getf(st, 'v');
        err = local_getf(st, 've');
        tdn = local_getf(st, 't');
end

if isempty(vel) || isempty(tdn), return; end

% Normalise orientation so rows = flowline points, cols = time.
vel = local_orient(vel, npts, numel(tdn));
err = local_orient(err, npts, numel(tdn));
vx  = local_orient(vx,  npts, numel(tdn));
vy  = local_orient(vy,  npts, numel(tdn));
tdn = tdn(:);
ok  = true;
end

% ---------------------------------------------------------------------------
function v = local_getf(st, f)
v = [];
if isstruct(st) && isfield(st, f)
    v = st.(f);
end
end

% ---------------------------------------------------------------------------
function A = local_orient(A, npts, nt)
%LOCAL_ORIENT  Ensure A is (npts x nt); transpose if stored as (nt x npts).
if isempty(A), return; end
if size(A,1)==npts
    return
elseif size(A,2)==npts && size(A,1)==nt
    A = A.';
elseif size(A,1)==nt && size(A,2)~=npts
    % Ambiguous but time is along dim1 -> transpose to put points on rows.
    A = A.';
end
end
