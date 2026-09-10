function FL = moc_flowline_to_struct(x, y, varargin)
%MOC_FLOWLINE_TO_STRUCT  Package a polyline into the toolkit's flowline layout.
%
%   Produces the same (x, y, d) geometry fields MOC_LOAD_OBS_FLOWLINE returns,
%   so the result is accepted by MOC_MODEL_FLOWLINE directly and gains
%   .vel/.time from MOC_SAMPLE_OBS_ALONG_FLOWLINE.
%
%   Usage:
%      FL = moc_flowline_to_struct(x, y)
%      FL = moc_flowline_to_struct(x, y, 'name','a045_06', 'spacing',200)
%      FL = moc_flowline_to_struct(x, y, 'extra',struct('glacier','a045'))
%
%   Input:
%      x, y : (np x 1) vertices [m], ordered from the downstream end inland
%   Name/value options:
%      'name'    : (char)   flowline name (default 'flowline')
%      'spacing' : (double) resample to this even along-line spacing [m];
%                  [] or 0 (default) keeps the shapefile's own vertices
%      'extra'   : (struct) extra fields to copy onto the output (e.g.
%                  glacier, flowline, terminus_date)
%
%   Output (struct FL):
%      .name   (char)      flowline name
%      .x, .y  (np x 1)    coordinates [m]
%      .d      (np x 1)    along-line distance [m], 0 at the downstream end
%      .length (double)    total length [m]
%      .epsg   (double)    projection code
%      plus whatever 'extra' carried

p = inputParser;
p.addParameter('name', 'flowline', @(s) ischar(s)||isstring(s));
p.addParameter('spacing', [], @(v) isempty(v) || isscalar(v));
p.addParameter('extra', struct(), @isstruct);
p.parse(varargin{:});
opt = p.Results;

x = x(:); y = y(:);
if numel(x) < 2 || numel(y) ~= numel(x)
    error('moc_flowline_to_struct:tooShort', ...
        'Need at least two matching (x,y) vertices to build a flowline.');
end

d = [0; cumsum(hypot(diff(x), diff(y)))];

% Drop repeated vertices, which would break the interp1 below
keep = [true; diff(d) > 0];
x = x(keep); y = y(keep); d = d(keep);
if numel(x) < 2
    error('moc_flowline_to_struct:degenerate', 'Flowline collapses to a point.');
end

if ~isempty(opt.spacing) && opt.spacing > 0
    dq = (0:opt.spacing:d(end))';
    if numel(dq) < 2, dq = [0; d(end)]; end
    x = interp1(d, x, dq, 'linear');
    y = interp1(d, y, dq, 'linear');
    d = dq;
end

cfg = config_moc();
FL = struct();
FL.name   = char(opt.name);
FL.x      = x;
FL.y      = y;
FL.d      = d;
FL.length = d(end);
FL.epsg   = cfg.epsg;

fn = fieldnames(opt.extra);
for i = 1:numel(fn)
    FL.(fn{i}) = opt.extra.(fn{i});
end

end
