function [xc, yc] = moc_clip_to_roi(x, y, roi, varargin)
%MOC_CLIP_TO_ROI  Clip a polyline to an ROI polygon.
%
%   A concave ROI (or one the line leaves and re-enters) can cut a line into
%   several pieces; this resolves them to a single connected line.
%
%   Usage:
%      [xc, yc] = moc_clip_to_roi(x, y, roi)
%      [xc, yc] = moc_clip_to_roi(x, y, roi, 'keep','first')
%
%   Input:
%      x, y : (np x 1) polyline vertices [m]
%      roi  : region of interest in any form MOC_ROI_POLYGON accepts
%   Name/value options:
%      'keep' : (char) 'longest'(default) keep the longest connected piece |
%               'first' keep the piece nearest the line's start
%
%   Output:
%      xc, yc : (nc x 1) the clipped line, empty if nothing falls inside

p = inputParser;
p.addParameter('keep', 'longest', @(s) ischar(s)||isstring(s));
p.parse(varargin{:});
keepMode = lower(char(p.Results.keep));

x = x(:); y = y(:);
xc = []; yc = [];
if numel(x) < 2, return; end

[rx, ry] = moc_roi_polygon(roi);
ws = warning('off', 'MATLAB:polyshape:repairedBySimplify');
cleanup = onCleanup(@() warning(ws));
pg = polyshape(rx(1:end-1), ry(1:end-1));
inLine = intersect(pg, [x y]);
if isempty(inLine), return; end

% Split the NaN-separated result into connected pieces
nanRow = isnan(inLine(:,1)) | isnan(inLine(:,2));
edges  = [0; find(nanRow); size(inLine,1)+1];
best = []; bestScore = -Inf;
for k = 1:numel(edges)-1
    seg = inLine(edges(k)+1:edges(k+1)-1, :);
    if size(seg,1) < 2, continue; end
    switch keepMode
        case 'longest'
            score = sum(hypot(diff(seg(:,1)), diff(seg(:,2))));
        case 'first'
            % earliest along the original line wins -> negate the distance
            score = -hypot(seg(1,1) - x(1), seg(1,2) - y(1));
        otherwise
            error('moc_clip_to_roi:badKeep', 'keep must be ''longest'' or ''first''.');
    end
    if score > bestScore, bestScore = score; best = seg; end
end
if isempty(best), return; end
xc = best(:,1); yc = best(:,2);

end
