function [xv, yv] = moc_roi_polygon(roi)
%MOC_ROI_POLYGON  Coerce anything ROI-shaped into closed polygon vertices.
%
%   Accepts every form the toolkit passes around for a region of interest and
%   returns one closed polygon in projected coordinates (cfg.epsg). Shapefile
%   coordinates are assumed to be in the same projection (both archives are
%   EPSG:3413) — nothing is reprojected here.
%
%   Usage:
%      [xv, yv] = moc_roi_polygon(roi)              % ROI struct (moc_load_roi)
%      [xv, yv] = moc_roi_polygon('sverdrups_trunk')% saved ROI name or path
%      [xv, yv] = moc_roi_polygon([xmin xmax ymin ymax])   % bbox, as roi.bbox
%      [xv, yv] = moc_roi_polygon([x(:) y(:)])      % np x 2 vertex matrix
%
%   Input:
%      roi : ROI struct with .x/.y | ROI name or file path (char/string) |
%            1x4 bbox [xmin xmax ymin ymax] | np x 2 matrix of vertices
%   Output:
%      xv, yv : (nv x 1 double) closed polygon vertices [m]

if ischar(roi) || isstring(roi)
    roi = moc_load_roi(char(roi));
end

if isstruct(roi)
    if ~isfield(roi,'x') || ~isfield(roi,'y')
        error('moc_roi_polygon:badRoi', 'ROI struct must have .x and .y vertices.');
    end
    xv = roi.x(:); yv = roi.y(:);
elseif isnumeric(roi) && isvector(roi) && numel(roi) == 4
    xmin = roi(1); xmax = roi(2); ymin = roi(3); ymax = roi(4);
    xv = [xmin; xmax; xmax; xmin; xmin];
    yv = [ymin; ymin; ymax; ymax; ymin];
elseif isnumeric(roi) && size(roi,2) == 2 && size(roi,1) >= 3
    xv = roi(:,1); yv = roi(:,2);
else
    error('moc_roi_polygon:badRoi', ...
        'Cannot interpret input as an ROI (struct, name, 1x4 bbox, or np x 2 matrix).');
end

% Close the ring if it is left open
if xv(1) ~= xv(end) || yv(1) ~= yv(end)
    xv(end+1) = xv(1);
    yv(end+1) = yv(1);
end

end
