function [xu, yu, scut] = moc_clip_upstream(x, y, tx, ty, varargin)
%MOC_CLIP_UPSTREAM  Cut a flowline at a terminus trace and keep the upstream part.
%
%   The terminus acts as a splitting feature: everything seaward of the
%   furthest-inland crossing is discarded, so a wiggly terminus that crosses
%   the flowline more than once leaves no seaward remnant.
%
%   Usage:
%      [xu, yu] = moc_clip_upstream(x, y, tx, ty)
%      [xu, yu] = moc_clip_upstream(x, y, tx, ty, 'downstream_end','last')
%
%   Input:
%      x, y   : (np x 1) flowline vertices [m]
%      tx, ty : (nt x 1) terminus trace vertices [m]
%   Name/value options:
%      'downstream_end' : (char) which end of the flowline is seaward —
%                         'auto'(default: whichever endpoint lies closer to the
%                         terminus; Felikson flowlines start at the terminus, so
%                         this resolves to 'first') | 'first' | 'last'
%
%   Output:
%      xu, yu : (nu x 1) the upstream part, starting exactly ON the terminus
%               and running inland. Empty if the two lines never cross.
%      scut   : (double) along-flowline distance of the cut [m] (NaN if none)

p = inputParser;
p.addParameter('downstream_end', 'auto', @(s) ischar(s)||isstring(s));
p.parse(varargin{:});
dsEnd = lower(char(p.Results.downstream_end));

x = x(:); y = y(:); tx = tx(:); ty = ty(:);
xu = []; yu = []; scut = NaN;
if numel(x) < 2 || numel(tx) < 2, return; end

if strcmp(dsEnd, 'auto')
    dFirst = local_point_to_polyline(x(1),   y(1),   tx, ty);
    dLast  = local_point_to_polyline(x(end), y(end), tx, ty);
    if dFirst <= dLast, dsEnd = 'first'; else, dsEnd = 'last'; end
end
switch dsEnd
    case 'first'
        % already terminus-first
    case 'last'
        x = flipud(x); y = flipud(y);
    otherwise
        error('moc_clip_upstream:badEnd', ...
            'downstream_end must be ''first'', ''last'' or ''auto''.');
end

[xi, yi] = polyxpoly(x, y, tx, ty);
if isempty(xi), return; end

% Cumulative distance along the flowline, and the arclength of each crossing
d = [0; cumsum(hypot(diff(x), diff(y)))];
s = arrayfun(@(k) local_project_onto_line(xi(k), yi(k), x, y, d), 1:numel(xi));
[scut, kcut] = max(s);
if ~isfinite(scut) || scut >= d(end), scut = NaN; return; end

% Upstream part: the crossing point itself, then every vertex beyond it
keep = d > scut;
xu = [xi(kcut); x(keep)];
yu = [yi(kcut); y(keep)];

end

% ===========================================================================
function s = local_project_onto_line(px, py, x, y, d)
%LOCAL_PROJECT_ONTO_LINE  Arclength of the point on the polyline nearest (px,py).
ax = x(1:end-1); ay = y(1:end-1);
bx = x(2:end);   by = y(2:end);
vx = bx - ax;    vy = by - ay;
len2 = vx.^2 + vy.^2;
t = ((px - ax).*vx + (py - ay).*vy) ./ max(len2, eps);
t = min(max(t, 0), 1);
qx = ax + t.*vx;  qy = ay + t.*vy;
[~, k] = min(hypot(px - qx, py - qy));
s = d(k) + t(k) * sqrt(len2(k));
end

% ---------------------------------------------------------------------------
function dmin = local_point_to_polyline(px, py, x, y)
%LOCAL_POINT_TO_POLYLINE  Shortest distance from a point to a polyline.
ax = x(1:end-1); ay = y(1:end-1);
bx = x(2:end);   by = y(2:end);
vx = bx - ax;    vy = by - ay;
t = ((px - ax).*vx + (py - ay).*vy) ./ max(vx.^2 + vy.^2, eps);
t = min(max(t, 0), 1);
dmin = min(hypot(px - (ax + t.*vx), py - (ay + t.*vy)));
end
