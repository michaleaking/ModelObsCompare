function F = moc_model_flowline(M, flowline, varargin)
%MOC_MODEL_FLOWLINE  Sample a model field along a flowline for every time step.
%
%   Convenience wrapper that samples the model at the vertices of a flowline
%   (from an observed flowline struct, an ROI/.exp contour, or explicit x/y) and
%   returns a distance-by-time matrix ready to compare against observed profiles.
%
%   Usage:
%      F = moc_model_flowline(M, Fobs)              % Fobs from moc_load_obs_flowline
%      F = moc_model_flowline(M, roi)               % use ROI polygon as the line
%      F = moc_model_flowline(M, [x y])             % np x 2 coordinate matrix
%      F = moc_model_flowline(M, flowline, 'field','vel')
%
%   Input:
%      M        : model struct from moc_load_model
%      flowline : observed-flowline struct (.x,.y[,.d]) | ROI struct | np x 2 matrix
%   Name/value options:
%      'field' : (char) 'vel'(default) | 'vx' | 'vy' | 'surface'
%
%   Output (struct F):
%      .x, .y  (np x 1)      flowline coordinates [m]
%      .d      (np x 1)      along-flowline distance [m] (0 at first point)
%      .time   (1 x nt)      decimal years
%      .val    (np x nt)     modelled field along the flowline
%      .field  (char)        field name

p = inputParser;
p.addParameter('field','vel', @(s) ischar(s)||isstring(s));
p.parse(varargin{:});
field = char(p.Results.field);

% Resolve flowline coordinates + distance
[x, y, d] = local_line_coords(flowline);

P = moc_model_at_points(M, x, y, 'field', field);

F = struct();
F.x = x; F.y = y; F.d = d;
F.time = P.time;
F.val = P.val;
F.field = field;

end

% ===========================================================================
function [x, y, d] = local_line_coords(fl)
d = [];
if isnumeric(fl) && size(fl,2)==2
    x = fl(:,1); y = fl(:,2);
elseif isstruct(fl) && isfield(fl,'x') && isfield(fl,'y')
    x = fl.x(:); y = fl.y(:);
    if isfield(fl,'d') && ~isempty(fl.d), d = fl.d(:); end
else
    error('moc_model_flowline:badLine', ...
        'flowline must be an obs struct, an ROI struct, or an np x 2 matrix.');
end
if isempty(d)
    % cumulative Euclidean distance along the line
    dseg = hypot(diff(x), diff(y));
    d = [0; cumsum(dseg)];
end
end
