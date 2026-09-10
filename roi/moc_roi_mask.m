function in = moc_roi_mask(roi, x, y)
%MOC_ROI_MASK  Test which points fall inside an ROI polygon.
%
%   Usage:
%      in = moc_roi_mask(roi, x, y)          % x,y same size -> logical same size
%      in = moc_roi_mask(roi, X, Y)          % works on grids from meshgrid too
%
%   Input:
%      roi  : ROI struct with polygon vertices .x, .y
%      x, y : (numeric arrays, matching size) query coordinates [m]
%   Output:
%      in   : (logical array, size of x) true where (x,y) is inside the ROI

if ~isfield(roi,'x') || ~isfield(roi,'y')
    error('moc_roi_mask:badRoi', 'roi must have .x and .y polygon vertices.');
end
in = inpolygon(x, y, roi.x, roi.y);

end
