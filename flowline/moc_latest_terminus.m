function [t, idx] = moc_latest_terminus(T, varargin)
%MOC_LATEST_TERMINUS  Most recent terminus trace, optionally one crossing a line.
%
%   Usage:
%      t = moc_latest_terminus(T)                       % newest in the set
%      t = moc_latest_terminus(T, 'line',[x y])         % newest one that crosses
%      t = moc_latest_terminus(T, 'line',[x y], 'max_distance',2000)
%
%   Input:
%      T : struct array from MOC_LOAD_TERMINI (has .date, .x, .y)
%   Name/value options:
%      'line'         : (np x 2) flowline coordinates; only traces that cross it
%                       are considered, falling back to the whole set if none do
%                       (so a flowline that has retreated past every trace still
%                       gets an answer)
%      'max_distance' : (double) when falling back, ignore traces further than
%                       this [m] from the line
%
%   Output:
%      t   : the newest matching element of T ([] if there is nothing to return)
%      idx : its index into the ORIGINAL T ([] if none)

p = inputParser;
p.addParameter('line', [], @(v) isempty(v) || (isnumeric(v) && size(v,2)==2));
p.addParameter('max_distance', [], @(v) isempty(v) || isscalar(v));
p.parse(varargin{:});
opt = p.Results;

t = []; idx = [];
if isempty(T), return; end

cand = 1:numel(T);
if ~isempty(opt.line)
    lx = opt.line(:,1); ly = opt.line(:,2);
    crosses = arrayfun(@(k) ~isempty(polyxpoly(lx, ly, T(k).x, T(k).y)), cand);
    if any(crosses)
        cand = cand(crosses);
    elseif ~isempty(opt.max_distance)
        near = arrayfun(@(k) local_line_gap(lx, ly, T(k).x, T(k).y), cand) ...
               <= opt.max_distance;
        if ~any(near), return; end
        cand = cand(near);
    end
end

[~, k] = max([T(cand).date]);
idx = cand(k);
t = T(idx);

end

% ===========================================================================
function g = local_line_gap(ax, ay, bx, by)
%LOCAL_LINE_GAP  Cheap vertex-to-vertex distance between two polylines.
g = Inf;
for i = 1:numel(bx)
    g = min(g, min(hypot(ax - bx(i), ay - by(i))));
end
end
