function T = moc_load_termini(roi, varargin)
%MOC_LOAD_TERMINI  Read traced terminus positions, keep those in an ROI.
%
%   Reads the Black & Joughin terminus shapefile (glacier_termini_v01.0.shp:
%   one LineString per glacier per image, with a SourceDate and a Quality_Fl)
%   and returns the traces inside an ROI, sorted oldest to newest.
%
%   Usage:
%      T = moc_load_termini(roi)
%      T = moc_load_termini(roi, 'max_quality',[])            % keep flagged traces
%      T = moc_load_termini(roi, 'trange',{'2019-01-01',''})  % bound the dates
%
%   Input:
%      roi : region of interest in any form MOC_ROI_POLYGON accepts, or []
%   Name/value options:
%      'path'        : (char)   terminus shapefile (default moc_path('termini_shp'))
%      'max_quality' : (double) keep traces with Quality_Fl <= this (default 0;
%                               the product flags lower-confidence traces with
%                               non-zero values). [] keeps all.
%      'trange'      : 1x2 datetime, cellstr or string of date bounds; either
%                      end may be empty/NaT for open-ended.
%
%   Output (1 x n struct array T, sorted by date):
%      .date       (datetime) SourceDate
%      .quality    (double)   Quality_Fl
%      .glacier_id (double)   Glacier_ID  (NOTE: unrelated to the flowline
%                             collection's glacier ids — pair them spatially)
%      .image_id   (char)     Image_ID
%      .x, .y      (np x 1)   trace vertices [m]

p = inputParser;
p.addParameter('path', '', @(s) ischar(s)||isstring(s));
p.addParameter('max_quality', 0, @(v) isempty(v) || isscalar(v));
p.addParameter('trange', [], @(v) isempty(v) || numel(v)==2);
p.parse(varargin{:});
opt = p.Results;

cfg = config_moc();
tpath = char(opt.path);
if isempty(tpath), tpath = moc_path('termini_shp'); end
if ~isfile(tpath)
    error('moc_load_termini:notFound', 'Terminus shapefile not found: %s', tpath);
end

haveRoi = ~isempty(roi);
if haveRoi
    [rx, ry] = moc_roi_polygon(roi);
    S = shaperead(tpath, 'BoundingBox', [min(rx) min(ry); max(rx) max(ry)]);
else
    S = shaperead(tpath);
end

T = struct('date',{},'quality',{},'glacier_id',{},'image_id',{},'x',{},'y',{});
for k = 1:numel(S)
    q = 0;
    if isfield(S, 'Quality_Fl'), q = double(S(k).Quality_Fl); end
    if ~isempty(opt.max_quality) && q > opt.max_quality, continue; end

    dt = local_parse_date(S(k));
    if isnat(dt), continue; end

    [x, y] = local_first_part(S(k).X, S(k).Y);
    if numel(x) < 2, continue; end
    if haveRoi && ~local_line_hits_roi(x, y, rx, ry), continue; end

    gid = NaN; if isfield(S,'Glacier_ID'), gid = double(S(k).Glacier_ID); end
    iid = '';  if isfield(S,'Image_ID'),   iid = char(string(S(k).Image_ID)); end
    T(end+1) = struct('date', dt, 'quality', q, 'glacier_id', gid, ...
        'image_id', iid, 'x', x, 'y', y); %#ok<AGROW>
end

if isempty(T), return; end

% Date bounds
if ~isempty(opt.trange)
    [lo, hi] = local_bounds(opt.trange);
    d = [T.date];
    keep = true(size(d));
    if ~isnat(lo), keep = keep & (d >= lo); end
    if ~isnat(hi), keep = keep & (d <= hi); end
    T = T(keep);
end
if isempty(T), return; end

[~, ord] = sort([T.date]);
T = T(ord);

end

% ===========================================================================
function dt = local_parse_date(rec)
%LOCAL_PARSE_DATE  SourceDate ('yyyy-MM-dd') to datetime, NaT if unparseable.
dt = NaT;
if ~isfield(rec, 'SourceDate'), return; end
s = strtrim(char(string(rec.SourceDate)));
if isempty(s), return; end
try
    dt = datetime(s, 'InputFormat', 'yyyy-MM-dd');
catch
    try, dt = datetime(s); catch, dt = NaT; end
end
end

% ---------------------------------------------------------------------------
function [lo, hi] = local_bounds(tr)
%LOCAL_BOUNDS  Turn a 1x2 date bound (datetime/cellstr/string) into two datetimes.
lo = NaT; hi = NaT;
if isdatetime(tr)
    lo = tr(1); hi = tr(2);
    return;
end
if isstring(tr), tr = cellstr(tr); end
if iscell(tr)
    if ~isempty(tr{1}), lo = datetime(tr{1}); end
    if ~isempty(tr{2}), hi = datetime(tr{2}); end
end
end

% ---------------------------------------------------------------------------
function [x, y] = local_first_part(X, Y)
%LOCAL_FIRST_PART  Column vectors for the first part of a shaperead record.
X = X(:); Y = Y(:);
n = find(isnan(X) | isnan(Y), 1, 'first');
if ~isempty(n), X = X(1:n-1); Y = Y(1:n-1); end
x = X; y = Y;
end

% ---------------------------------------------------------------------------
function tf = local_line_hits_roi(x, y, rx, ry)
%LOCAL_LINE_HITS_ROI  True if the trace touches the ROI polygon at all.
tf = any(inpolygon(x, y, rx, ry));
if ~tf
    tf = ~isempty(polyxpoly(x, y, rx, ry));
end
end
