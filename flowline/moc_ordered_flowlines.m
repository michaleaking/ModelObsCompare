function [idx, F] = moc_ordered_flowlines(F, varargin)
%MOC_ORDERED_FLOWLINES  Order a glacier's fan of flowlines from the centre out.
%
%   Felikson's flowlines are seeded side by side across the trunk, so ranking
%   them by how far their downstream end sits from the mean of all the
%   downstream ends walks the fan centre first, then alternately out to either
%   side. That order is what MOC_FLOWLINE_FROM_ROI follows when the centre line
%   has no terminus trace crossing it.
%
%   Usage:
%      idx = moc_ordered_flowlines(F)                    % over the pooled set
%      idx = moc_ordered_flowlines(F, 'glacier','a045')  % one glacier's fan
%
%   Input:
%      F : struct array from MOC_LOAD_FLOWLINES
%   Name/value options:
%      'glacier' : (char) restrict to this glacier id (e.g. 'a045')
%
%   Output:
%      idx : (1 x n) indices into the ORIGINAL F, centre first
%      F   : the input, returned for convenience (F(idx) is the ordered fan)

p = inputParser;
p.addParameter('glacier', '', @(s) ischar(s)||isstring(s));
p.parse(varargin{:});
gid = char(p.Results.glacier);

if isempty(F)
    error('moc_ordered_flowlines:empty', 'No flowlines to choose from.');
end

sel = true(1, numel(F));
if ~isempty(gid)
    sel = strcmp({F.glacier}, gid);
    if ~any(sel)
        error('moc_ordered_flowlines:noGlacier', ...
            'Glacier "%s" has no flowlines here (available: %s).', ...
            gid, strjoin(unique({F.glacier}), ', '));
    end
end
which = find(sel);

if numel(which) == 1
    idx = which; return;
end

% Downstream ends (point 1 of each line), ranked by distance from their mean
ends = [arrayfun(@(k) F(k).x(1), which)', arrayfun(@(k) F(k).y(1), which)'];
ctr  = mean(ends, 1);
[~, ord] = sort(hypot(ends(:,1) - ctr(1), ends(:,2) - ctr(2)));
idx = which(ord);

end
