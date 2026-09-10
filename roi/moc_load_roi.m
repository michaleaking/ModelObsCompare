function roi = moc_load_roi(name, varargin)
%MOC_LOAD_ROI  Load a previously saved ROI (.mat) or build one from an .exp contour.
%
%   Usage:
%      roi = moc_load_roi('upernavik')            % rois/upernavik.mat or .exp
%      roi = moc_load_roi('/path/to/contour.exp')
%
%   Input:
%      name : (char) ROI name (looked up in config_moc().roi_dir) or a full path
%             to a .mat or .exp file.
%   Name/value options:
%      'dir' : (char) folder to search (default: config_moc().roi_dir)
%
%   Output:
%      roi : ROI struct with .name, .type, .x, .y, .bbox, .epsg

p = inputParser;
p.addParameter('dir','', @(s) ischar(s)||isstring(s));
p.parse(varargin{:});
odir = char(p.Results.dir);
if isempty(odir)
    odir = moc_path('roi_dir');
end
name = char(name);

% Resolve to a concrete file
candidates = {name, ...
    fullfile(odir, name), ...
    fullfile(odir, [name '.mat']), ...
    fullfile(odir, [name '.exp'])};
target = '';
for i = 1:numel(candidates)
    if isfile(candidates{i}), target = candidates{i}; break; end
end
if isempty(target)
    error('moc_load_roi:notFound', 'No ROI found for "%s" in %s', name, odir);
end

[~, base, ext] = fileparts(target);
switch lower(ext)
    case '.mat'
        S = load(target);
        roi = S.roi;
    case '.exp'
        c = expread(target);
        roi = struct();
        roi.name = base;
        roi.type = 'polygon';
        roi.x = c(1).x(:);
        roi.y = c(1).y(:);
        roi.bbox = [min(roi.x) max(roi.x) min(roi.y) max(roi.y)];
        cfg = config_moc();
        roi.epsg = cfg.epsg;
        roi.file = target;
    otherwise
        error('moc_load_roi:badExt', 'Unsupported ROI file type: %s', ext);
end

end
