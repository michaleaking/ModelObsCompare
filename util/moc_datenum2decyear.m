function yr = moc_datenum2decyear(dn)
%MOC_DATENUM2DECYEAR  Convert MATLAB datenum to decimal year.
%
%   Observed flowline files (redrawn_Jan25) store time as MATLAB datenum,
%   whereas ISSM model output uses decimal years. Use this to align them.
%
%   Usage:
%      yr = moc_datenum2decyear(dn)
%
%   Input:
%      dn : (numeric array) MATLAB serial date numbers
%   Output:
%      yr : (double array, same shape) decimal years (e.g. 2018.5)

sz = size(dn);
dn = dn(:);
dv = datevec(dn);
y  = dv(:,1);
yearStart = datenum(y,     1, 1);
yearEnd   = datenum(y + 1, 1, 1);
yr = y + (dn - yearStart) ./ (yearEnd - yearStart);
yr = reshape(yr, sz);

end
