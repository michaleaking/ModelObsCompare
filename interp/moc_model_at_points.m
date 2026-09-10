function P = moc_model_at_points(M, x, y, varargin)
%MOC_MODEL_AT_POINTS  Interpolate a model field onto (x,y) points for every time step.
%
%   Uses ISSM's InterpFromMesh2d to sample the model mesh at arbitrary points.
%   Returns a points-by-time matrix, giving a modelled time series at each point
%   (e.g. along a flowline or at a single station coordinate).
%
%   Usage:
%      P = moc_model_at_points(M, x, y)
%      P = moc_model_at_points(M, x, y, 'field','vel')
%
%   Input:
%      M    : model struct from moc_load_model
%      x, y : (np x 1) query coordinates [m] (same EPSG as the model mesh)
%   Name/value options:
%      'field' : (char) 'vel'(default) | 'vx' | 'vy' | 'surface'
%
%   Output (struct P):
%      .x, .y  (np x 1)         query points
%      .time   (1 x nt)         decimal years (NaN for a snapshot)
%      .val    (np x nt)        interpolated field at each point/time
%      .field  (char)           field name

p = inputParser;
p.addParameter('field','vel', @(s) ischar(s)||isstring(s));
p.parse(varargin{:});
field = char(p.Results.field);

if ~isfield(M, field)
    error('moc_model_at_points:badField', 'Model struct has no field "%s".', field);
end
x = x(:); y = y(:);
np = numel(x);
nt = size(M.(field), 2);

val = nan(np, nt);
for t = 1:nt
    val(:, t) = InterpFromMesh2d(M.elements, M.x, M.y, M.(field)(:, t), x, y);
end

P = struct();
P.x = x; P.y = y;
P.time = M.time;
P.val = val;
P.field = field;

end
