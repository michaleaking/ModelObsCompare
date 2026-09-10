function O = moc_load_obs_netcdf(ncpath, varargin)
%MOC_LOAD_OBS_NETCDF  Load a gridded observed-velocity netCDF (Sentinel/CSK subset).
%
%   These files store a 4-D `velocity(time, component, y, x)` array where the
%   component dimension holds the bands {'vv','vx','vy','ex','ey','dT'}. This
%   loader is robust to dimension ORDER: it locates each axis by matching sizes
%   against the x/y/time/component coordinate variables, then permutes each band
%   to a consistent (ny x nx x nt) layout.
%
%   Usage:
%      O = moc_load_obs_netcdf(ncpath)
%      O = moc_load_obs_netcdf('Sentinel_Subset_Upernavik.nc')   % bare name -> obs_netcdf_dir
%      O = moc_load_obs_netcdf(ncpath, 'trange',[2018 2020], 'bands',{'vx','vy','vv'})
%
%   Input:
%      ncpath : (char) path or bare filename (resolved against obs_netcdf_dir)
%   Name/value options:
%      'varname' : (char)    velocity variable name (default 'velocity')
%      'bands'   : (cellstr) subset of bands to return (default: all present)
%      'trange'  : (1x2)     [t0 t1] decimal-year window to keep (default: all)
%
%   Output (struct O):
%      .file          (char)              source file
%      .x, .y         (1xnx / 1xny)       projected coordinates [m]
%      .time          (1 x nt double)     decimal years
%      .bands         (cellstr)           band names returned
%      .(band)        (ny x nx x nt)      one field per band (e.g. O.vx, O.vy, O.vv)
%      .epsg          (double)            projection code parsed from spatial_ref (best effort)
%
%   Nodata handling: these products declare a useless `_FillValue` ATTRIBUTE
%   (NaN) and keep the real per-band sentinels in a `_FillValue` VARIABLE, e.g.
%   {-1, -2e9, -2e9, -1, -1} for {vv, vx, vy, ex, ey}. Those sentinels are set
%   to NaN here, along with any non-physical magnitude (>=1e9) and, on the
%   physically non-negative bands (vv/ex/ey/dT), any negative value. Without
%   this a nodata pixel reads as a -1 m/yr "speed". Same rules as the Python
%   reader (moc_py.obs_netcdf).

p = inputParser;
p.addParameter('varname', 'velocity', @(s) ischar(s)||isstring(s));
p.addParameter('bands', {}, @(c) iscell(c) || isstring(c));
p.addParameter('trange', [], @(v) isempty(v) || (isnumeric(v) && numel(v)==2));
p.parse(varargin{:});
opt = p.Results;
varname = char(opt.varname);

% Resolve path
if ~isfile(ncpath)
    cfg = config_moc();
    cand = fullfile(moc_path('obs_netcdf_dir'), ncpath);
    if isfile(cand), ncpath = cand; else
        error('moc_load_obs_netcdf:notFound', 'netCDF not found: %s', ncpath);
    end
end

info = ncinfo(ncpath);

% --- coordinate vectors ---
x = double(ncread(ncpath, 'x'));  x = x(:)';
y = double(ncread(ncpath, 'y'));  y = y(:)';
tRaw = double(ncread(ncpath, 'time'));  tRaw = tRaw(:)';
tUnits = local_get_att(info, 'time', 'units');   % e.g. 'days since 2015-01-06 12:00:00'
tdn = local_time_to_datenum(tRaw, tUnits);
tYear = moc_datenum2decyear(tdn);

% --- component / band names ---
bandNames = local_component_names(ncpath, info);
nComp = numel(bandNames);

% --- read velocity and figure out its dimension order ---
vinfo = local_varinfo(info, varname);
vsize = vinfo.Size;
% Match each dimension of `velocity` to x/y/time/component by size.
role = local_map_dims(vsize, numel(x), numel(y), numel(tRaw), nComp);

% Which bands to return
wantBands = opt.bands;
if isempty(wantBands), wantBands = bandNames; end
wantBands = cellstr(wantBands);

% Optional time subset -> use ncread start/count on the time dimension only
tKeep = true(1, numel(tYear));
if ~isempty(opt.trange)
    tKeep = tYear >= opt.trange(1) & tYear <= opt.trange(2);
    if ~any(tKeep)
        warning('moc_load_obs_netcdf:noTimes', 'No frames within trange; returning empty.');
    end
end

O = struct();
O.file  = ncpath;
O.x     = x;
O.y     = y;
O.time  = tYear(tKeep);
O.bands = wantBands;
O.epsg  = local_parse_epsg(info);

% Per-band nodata sentinels, if the file carries them (see the note above)
fillVals = local_fill_values(ncpath, nComp);

% Read each requested band, permute to (ny, nx, nt)
for b = 1:numel(wantBands)
    ci = find(strcmp(bandNames, wantBands{b}), 1);
    if isempty(ci)
        warning('moc_load_obs_netcdf:badBand', 'Band "%s" not present; skipping.', wantBands{b});
        continue
    end
    % start/count over the raw variable dimensions
    start = ones(1, numel(vsize));
    count = vsize;                 % read all, then slice component (avoids per-frame reads)
    start(role.comp) = ci;
    count(role.comp) = 1;
    raw = ncread(ncpath, varname, start, count);   % singleton on comp dim
    raw = squeeze(raw);
    % After squeeze, remaining dims are {x,y,t} in their original relative order.
    fld = local_to_yxt(raw, vsize, role);
    fld = fld(:, :, tKeep);
    O.(wantBands{b}) = local_mask_fill(double(fld), fillVals(ci), wantBands{b});
end

fprintf('Loaded obs grid: %s | %dx%d px, %d frames (%.2f-%.2f), bands: %s\n', ...
    O.file, numel(x), numel(y), numel(O.time), min(O.time), max(O.time), strjoin(wantBands, ','));

end

% ===========================================================================
function names = local_component_names(ncpath, info)
%LOCAL_COMPONENT_NAMES  Read the 'component' band labels, falling back to config order.
names = {};
try
    c = ncread(ncpath, 'component');
    if iscell(c)
        names = cellfun(@char, c, 'UniformOutput', false);
    elseif ischar(c)
        names = cellstr(c);
    elseif isstring(c)
        names = cellstr(c);
    end
catch
end
names = names(:)';
if isempty(names)
    cfg = config_moc();
    names = cfg.obs_bands;   % documented default order
end
end

% ---------------------------------------------------------------------------
function fv = local_fill_values(ncpath, nComp)
%LOCAL_FILL_VALUES  Per-band nodata sentinels from the '_FillValue' VARIABLE.
%   Returns a 1 x nComp vector of sentinels (NaN where none is declared). The
%   like-named attribute on `velocity` is NaN in these files and is no help.
fv = NaN(1, nComp);
try
    v = double(ncread(ncpath, '_FillValue'));
    v = v(:)';
    n = min(numel(v), nComp);
    fv(1:n) = v(1:n);
catch
    % no such variable — the rule-based masking below still applies
end
end

% ---------------------------------------------------------------------------
function fld = local_mask_fill(fld, fillval, band)
%LOCAL_MASK_FILL  Set nodata sentinels to NaN.
%   Catches the declared per-band sentinel, any non-physical magnitude (>=1e9,
%   e.g. the -2e9 used on vx/vy), and negatives on the bands that cannot
%   physically be negative (vv is a speed magnitude, ex/ey are error
%   magnitudes, dT an image-pair separation — these files flag those with -1).
if ~isempty(fillval) && isfinite(fillval)
    fld(fld == fillval) = NaN;
end
fld(abs(fld) >= 1e9) = NaN;
if any(strcmp(band, {'vv','ex','ey','dT'}))
    fld(fld < 0) = NaN;
end
end

% ---------------------------------------------------------------------------
function vinfo = local_varinfo(info, varname)
idx = find(strcmp({info.Variables.Name}, varname), 1);
if isempty(idx)
    error('moc_load_obs_netcdf:noVar', 'Variable "%s" not in file.', varname);
end
vinfo = info.Variables(idx);
end

% ---------------------------------------------------------------------------
function role = local_map_dims(vsize, nx, ny, nt, nc)
%LOCAL_MAP_DIMS  Identify which velocity dimension is x / y / time / component by size.
role.x = find(vsize==nx, 1);
role.y = find(vsize==ny, 1);
role.t = find(vsize==nt, 1);
role.comp = find(vsize==nc, 1);
if any(cellfun(@isempty, {role.x, role.y, role.t, role.comp}))
    error('moc_load_obs_netcdf:dimMatch', ...
        ['Could not unambiguously match velocity dims %s to x=%d y=%d t=%d comp=%d. ' ...
         'Sizes may be non-unique; specify a different varname or inspect with ncinfo.'], ...
        mat2str(vsize), nx, ny, nt, nc);
end
end

% ---------------------------------------------------------------------------
function fld = local_to_yxt(raw, vsize, role)
%LOCAL_TO_YXT  Permute a single-band array (component squeezed out) to (ny, nx, nt).
% After squeezing out the component dim, the surviving dims keep their original
% relative order; find where y, x, t landed and permute to (y, x, t).
dims = setdiff(1:numel(vsize), role.comp, 'stable');   % raw dims after removing comp
% position of y,x,t within the squeezed array
py = find(dims==role.y);
px = find(dims==role.x);
pt = find(dims==role.t);
fld = permute(raw, [py, px, pt]);
end

% ---------------------------------------------------------------------------
function s = local_get_att(info, varName, attName)
s = '';
vi = find(strcmp({info.Variables.Name}, varName), 1);
if isempty(vi), return; end
atts = info.Variables(vi).Attributes;
if isempty(atts), return; end
ai = find(strcmp({atts.Name}, attName), 1);
if ~isempty(ai), s = atts(ai).Value; end
end

% ---------------------------------------------------------------------------
function dn = local_time_to_datenum(traw, units)
%LOCAL_TIME_TO_DATENUM  Convert CF "<step> since <ref>" times to MATLAB datenum.
dn = traw;   % fallback: assume already datenum-like
if isempty(units), return; end
tok = regexp(units, '(\w+)\s+since\s+(.+)', 'tokens', 'once');
if isempty(tok), return; end
step = lower(strtrim(tok{1}));
refStr = strtrim(tok{2});
refStr = strrep(refStr, 'T', ' ');           % ISO 'T' separator -> space
refStr = regexprep(refStr, '\s*(UTC|Z)$', ''); % drop trailing tz marker
ref = NaN;
fmts = {'yyyy-mm-dd HH:MM:SS','yyyy-mm-dd HH:MM','yyyy-mm-dd'};
for f = 1:numel(fmts)
    try, ref = datenum(refStr, fmts{f}); break; catch, end
end
if isnan(ref), ref = datenum(refStr); end     % last resort
switch step
    case {'days','day'},        dn = ref + traw;
    case {'hours','hour'},      dn = ref + traw/24;
    case {'minutes','minute'},  dn = ref + traw/1440;
    case {'seconds','second'},  dn = ref + traw/86400;
    otherwise,                  dn = ref + traw;  % assume days
end
end

% ---------------------------------------------------------------------------
function epsg = local_parse_epsg(info)
epsg = NaN;
vi = find(strcmp({info.Variables.Name}, 'spatial_ref'), 1);
if isempty(vi), return; end
atts = info.Variables(vi).Attributes;
if isempty(atts), return; end
ai = find(strcmp({atts.Name}, 'spatial_ref'), 1);
if isempty(ai), ai = find(strcmp({atts.Name}, 'crs_wkt'), 1); end
if isempty(ai), return; end
tok = regexp(atts(ai).Value, 'EPSG","(\d+)"\]\]$', 'tokens', 'once');
if ~isempty(tok), epsg = str2double(tok{1}); end
end
