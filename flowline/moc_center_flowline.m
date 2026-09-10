function [c, idx] = moc_center_flowline(F, varargin)
%MOC_CENTER_FLOWLINE  Pick the middle flowline of a glacier's fan of flowlines.
%
%   The "centre" line is the one whose downstream end sits closest to the mean
%   of all the downstream ends — the first of MOC_ORDERED_FLOWLINES.
%
%   Usage:
%      c = moc_center_flowline(F)                    % centre of the pooled set
%      c = moc_center_flowline(F, 'glacier','a045')  % centre of one glacier's fan
%      [c, idx] = moc_center_flowline(F, ...)        % also its index into F
%
%   Input:
%      F : struct array from MOC_LOAD_FLOWLINES
%   Name/value options:
%      'glacier' : (char) restrict to this glacier id (e.g. 'a045')
%
%   Output:
%      c   : the selected element of F
%      idx : its index into the ORIGINAL F
%
%   See also MOC_ORDERED_FLOWLINES.

p = inputParser;
p.addParameter('glacier', '', @(s) ischar(s)||isstring(s));
p.parse(varargin{:});
gid = char(p.Results.glacier);

ord = moc_ordered_flowlines(F, 'glacier', gid);
idx = ord(1);
c = F(idx);

end
