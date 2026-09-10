function FL = moc_flowline_from_roi(roi, varargin)
%MOC_FLOWLINE_FROM_ROI  One ROI-cropped, terminus-clipped flowline, ready to compare.
%
%   Runs the whole pipeline: read every flowline shapefile, keep those reaching
%   into the ROI, take the centre flowline of the glacier, split it at the most
%   recent terminus trace that crosses it, keep the upstream side, clip that to
%   the ROI, and package it.
%
%   If no terminus trace crosses the centre flowline, the next line out is
%   tried instead — working outwards from the centre, alternating sides — and
%   the first one a trace does cross is used. Compare the .flowline and
%   .center_flowline fields to see whether that happened. Naming a 'flowline'
%   explicitly disables the search: that line is used or the call fails.
%
%   Usage:
%      FL = moc_flowline_from_roi(roi)
%      FL = moc_flowline_from_roi('sverdrups_trunk', 'spacing',200)
%      FL = moc_flowline_from_roi(roi, 'glacier','b045', 'flowline','05')
%
%   Input:
%      roi : region of interest in any form MOC_ROI_POLYGON accepts
%   Name/value options:
%      'glacier'          : (char)   glacier id (default: the glacier with the
%                           most flowlines in the ROI)
%      'flowline'         : (char)   flowline id (default: the centre one,
%                           falling outwards to its neighbours if no terminus
%                           trace crosses it)
%      'spacing'          : (double) even resampling interval [m] ([] = native)
%      'name'             : (char)   name for the output ('<glacier>_<flowline>')
%      'dir'              : (char)   flowline shapefile folder
%      'termini_path'     : (char)   terminus shapefile
%      'max_quality'      : (double) terminus quality cut-off (default 0)
%      'trange'           : 1x2 date bounds on the terminus traces — use it to
%                           clip at the front position of a particular epoch
%      'require_terminus' : (logical) error if no trace crosses ANY candidate
%                           flowline (default true; false falls back to the
%                           un-split centre line)
%
%   Output (struct FL): as MOC_FLOWLINE_TO_STRUCT, plus
%      .glacier, .flowline (the line actually used), .center_flowline (the
%      centre of the fan, which differs when the search fell outwards),
%      .terminus_date (char, '' if not clipped), .terminus_file, .file
%
%   See also MOC_ORDERED_FLOWLINES, MOC_SAMPLE_OBS_ALONG_FLOWLINE,
%   MOC_MODEL_FLOWLINE.

p = inputParser;
p.addParameter('glacier', '', @(s) ischar(s)||isstring(s));
p.addParameter('flowline', '', @(s) ischar(s)||isstring(s));
p.addParameter('spacing', [], @(v) isempty(v) || isscalar(v));
p.addParameter('name', '', @(s) ischar(s)||isstring(s));
p.addParameter('dir', '', @(s) ischar(s)||isstring(s));
p.addParameter('termini_path', '', @(s) ischar(s)||isstring(s));
p.addParameter('max_quality', 0, @(v) isempty(v) || isscalar(v));
p.addParameter('trange', [], @(v) isempty(v) || numel(v)==2);
p.addParameter('require_terminus', true, @(v) islogical(v)||isnumeric(v));
p.parse(varargin{:});
opt = p.Results;

cfg = config_moc();
tpath = char(opt.termini_path);
if isempty(tpath), tpath = moc_path('termini_shp'); end

% --- flowlines in the ROI ------------------------------------------------
F = moc_load_flowlines(roi, 'dir', char(opt.dir));
if isempty(F)
    error('moc_flowline_from_roi:noFlowlines', 'No flowlines intersect this ROI.');
end

gid = char(opt.glacier);
if isempty(gid)
    g = {F.glacier};
    u = unique(g);
    counts = cellfun(@(a) sum(strcmp(g, a)), u);
    [~, k] = max(counts);
    gid = u{k};
end

ranked = moc_ordered_flowlines(F, 'glacier', gid);   % centre first, then outwards
centerId = F(ranked(1)).flowline;

flid = char(opt.flowline);
if isempty(flid)
    cand = ranked;
else
    sel = strcmp({F.glacier}, gid) & strcmp({F.flowline}, flid);
    if ~any(sel)
        avail = {F(strcmp({F.glacier}, gid)).flowline};
        error('moc_flowline_from_roi:noFlowline', ...
            'Flowline "%s" not in glacier "%s" (available: %s).', ...
            flid, gid, strjoin(avail, ', '));
    end
    cand = find(sel, 1, 'first');       % named explicitly: no outward search
end

T = moc_load_termini(roi, 'path', tpath, ...
    'max_quality', opt.max_quality, 'trange', opt.trange);

% --- walk out from the centre until a trace actually crosses a line ------
row = []; x = []; y = []; tdate = '';
for k = cand
    term = moc_latest_terminus(T, 'line', [F(k).x F(k).y]);
    if isempty(term), continue; end
    [xu, yu] = moc_clip_upstream(F(k).x, F(k).y, term.x, term.y);
    if isempty(xu), continue; end
    row = F(k); x = xu; y = yu;
    tdate = datestr(term.date, 'yyyy-mm-dd'); %#ok<DATST>
    break;
end

if isempty(row)
    if opt.require_terminus
        error('moc_flowline_from_roi:noTerminus', ...
            ['No terminus trace in the ROI crosses any flowline of glacier ' ...
             '"%s" (tried %s); pass ''require_terminus'',false to skip the clip.'], ...
            gid, strjoin({F(cand).flowline}, ', '));
    end
    row = F(cand(1));                  % un-split, from the centre outwards
    x = row.x; y = row.y;
end

% --- clip to the ROI and package -----------------------------------------
[x, y] = moc_clip_to_roi(x, y, roi);
if isempty(x)
    error('moc_flowline_from_roi:emptyAfterClip', ...
        ['Nothing of the flowline is left inside the ROI once the part ' ...
         'downstream of the terminus is removed.']);
end

name = char(opt.name);
if isempty(name), name = sprintf('%s_%s', gid, row.flowline); end

extra = struct('glacier', gid, 'flowline', row.flowline, ...
    'center_flowline', centerId, ...
    'terminus_date', tdate, 'terminus_file', tpath, 'file', row.file);
FL = moc_flowline_to_struct(x, y, 'name', name, ...
    'spacing', opt.spacing, 'extra', extra);

end
