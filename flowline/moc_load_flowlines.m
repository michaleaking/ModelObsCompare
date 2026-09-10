function F = moc_load_flowlines(roi, varargin)
%MOC_LOAD_FLOWLINES  Read every per-glacier flowline shapefile, keep those in an ROI.
%
%   Loops the whole Felikson et al. (2020) collection (glacier0001.shp ...
%   glacierd216.shp), each file holding several flowlines identified by a
%   `flowline` attribute, and returns them as one struct array. Flowline
%   vertices run from the terminus inland, so point 1 is the downstream end.
%
%   Usage:
%      F = moc_load_flowlines(roi)                      % ROI struct/name/bbox
%      F = moc_load_flowlines([])                       % the whole collection
%      F = moc_load_flowlines(roi, 'clip',true)         % trim geometry to the ROI
%      F = moc_load_flowlines(roi, 'iterations',true)   % keep *_iterNN.shp too
%
%   Input:
%      roi : region of interest in any form MOC_ROI_POLYGON accepts, or []
%            for no spatial filter
%   Name/value options:
%      'dir'        : (char)    folder of shapefiles (default moc_path('flowline_shp_dir'))
%      'pattern'    : (char)    glob for the files (default 'glacier*.shp')
%      'iterations' : (logical) keep the *_iterNN.shp intermediate versions (false)
%      'clip'       : (logical) clip geometry to the ROI instead of returning
%                               whole flowlines that merely reach into it (false)
%      'predicate'  : (char)    'intersects'(default) any overlap | 'within'
%                               entirely inside the ROI
%
%   Output (1 x n struct array F):
%      .glacier   (char)   glacier id parsed from the filename, e.g. 'a045'
%      .flowline  (char)   flowline id within that glacier, e.g. '05'
%      .iteration (double) 0 for the final flowline, N for *_iterNN
%      .file      (char)   source shapefile
%      .x, .y     (np x 1) vertices [m], terminus end first

p = inputParser;
p.addParameter('dir', '', @(s) ischar(s)||isstring(s));
p.addParameter('pattern', 'glacier*.shp', @(s) ischar(s)||isstring(s));
p.addParameter('iterations', false, @(v) islogical(v)||isnumeric(v));
p.addParameter('clip', false, @(v) islogical(v)||isnumeric(v));
p.addParameter('predicate', 'intersects', @(s) ischar(s)||isstring(s));
p.parse(varargin{:});
opt = p.Results;
predicate = lower(char(opt.predicate));

cfg = config_moc();
fdir = char(opt.dir);
if isempty(fdir), fdir = moc_path('flowline_shp_dir'); end

files = dir(fullfile(fdir, char(opt.pattern)));
if isempty(files)
    error('moc_load_flowlines:noFiles', ...
        'No shapefiles matching %s in %s', char(opt.pattern), fdir);
end

haveRoi = ~isempty(roi);
if haveRoi
    [rx, ry] = moc_roi_polygon(roi);
    bb = [min(rx) min(ry); max(rx) max(ry)];   % shaperead prefilter
end

F = struct('glacier',{},'flowline',{},'iteration',{},'file',{},'x',{},'y',{});
for i = 1:numel(files)
    [~, base] = fileparts(files(i).name);
    [gid, iter, isFlowFile] = local_parse_name(base);
    if ~isFlowFile, continue; end                 % e.g. sverdrup_masks.shp
    if iter > 0 && ~opt.iterations, continue; end

    fpath = fullfile(fdir, files(i).name);
    if haveRoi
        S = shaperead(fpath, 'BoundingBox', bb);
    else
        S = shaperead(fpath);
    end
    if isempty(S) || ~isfield(S, 'flowline'), continue; end

    for k = 1:numel(S)
        [x, y] = local_first_part(S(k).X, S(k).Y);
        if numel(x) < 2, continue; end
        if haveRoi
            if ~local_line_in_roi(x, y, rx, ry, predicate), continue; end
            if opt.clip
                [x, y] = moc_clip_to_roi(x, y, [rx ry]);
                if isempty(x), continue; end
            end
        end
        F(end+1) = struct('glacier', gid, 'flowline', char(string(S(k).flowline)), ...
            'iteration', iter, 'file', fpath, 'x', x, 'y', y); %#ok<AGROW>
    end
end

end

% ===========================================================================
function [gid, iter, ok] = local_parse_name(base)
%LOCAL_PARSE_NAME  Pull the glacier id and iteration number out of a filename.
%   'glaciera045'        -> gid 'a045', iter 0
%   'glaciera045_iter01' -> gid 'a045', iter 1
%   anything else        -> ok = false  (e.g. sverdrup_masks.shp)
%   Parsed by hand rather than with one regexp: MATLAB's regexp drops the
%   token of an optional group, so '_iter01' would silently read as iter 0.
gid = ''; iter = 0; ok = false;
prefix = 'glacier';
if ~startsWith(base, prefix), return; end
rest = base(numel(prefix)+1:end);

k = strfind(rest, '_iter');
if ~isempty(k)
    iter = str2double(rest(k(1)+5:end));
    if isnan(iter), iter = 0; return; end
    rest = rest(1:k(1)-1);
end
if isempty(regexp(rest, '^[a-z]?\d{3,4}$', 'once')), iter = 0; return; end
gid = rest;
ok = true;
end

% ---------------------------------------------------------------------------
function [x, y] = local_first_part(X, Y)
%LOCAL_FIRST_PART  Column vectors for the first part of a shaperead record.
%   shaperead terminates each feature (and separates multi-part features) with
%   NaN; keeping only the first part avoids joining disjoint pieces into one line.
X = X(:); Y = Y(:);
n = find(isnan(X) | isnan(Y), 1, 'first');
if ~isempty(n), X = X(1:n-1); Y = Y(1:n-1); end
x = X; y = Y;
end

% ---------------------------------------------------------------------------
function tf = local_line_in_roi(x, y, rx, ry, predicate)
%LOCAL_LINE_IN_ROI  Spatial test of a polyline against the ROI polygon.
in = inpolygon(x, y, rx, ry);
switch predicate
    case 'within'
        tf = all(in);
    case 'intersects'
        tf = any(in);
        if ~tf   % no vertex inside, but the line may still cut across
            xi = polyxpoly(x, y, rx, ry);
            tf = ~isempty(xi);
        end
    otherwise
        error('moc_load_flowlines:badPredicate', ...
            'predicate must be ''intersects'' or ''within''.');
end
end
