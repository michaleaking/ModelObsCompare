function T = moc_list_models(varargin)
%MOC_LIST_MODELS  List available ISSM model output files with parsed metadata.
%
%   Scans config_moc().model_dir for .mat files and parses common tokens out
%   of the ISSM naming convention (region, fric law, catchment, years) so you
%   can pick a run without remembering exact filenames.
%
%   Usage:
%      T = moc_list_models()
%      T = moc_list_models('dir', '/some/other/folder')
%
%   Name/value options:
%      'dir' : (char) folder to scan (default: config_moc().model_dir)
%
%   Output:
%      T : table with columns file, name, region, fric, catchment, bytes

p = inputParser;
p.addParameter('dir', '', @(s) ischar(s) || isstring(s));
p.parse(varargin{:});
mdir = char(p.Results.dir);
if isempty(mdir)
    cfg = config_moc();
    mdir = moc_path('model_dir');
end

d = dir(fullfile(mdir, 'Model_*.mat'));
% Some model files were saved without a .mat extension; include them too.
d = [d; dir(fullfile(mdir, 'Model_*'))];
[~, iu] = unique({d.name});
d = d(iu);
d = d(~[d.isdir]);

if isempty(d)
    warning('moc_list_models:empty', 'No Model_* files found in %s', mdir);
    T = table();
    return
end

n = numel(d);
name = strings(n,1); region = strings(n,1);
fric = nan(n,1); catchment = nan(n,1); bytes = nan(n,1);

for i = 1:n
    name(i)  = string(d(i).name);
    bytes(i) = d(i).bytes;
    tok = regexp(d(i).name, 'Model_(NW|CW)', 'tokens', 'once');
    if ~isempty(tok), region(i) = string(tok{1}); end
    tok = regexp(d(i).name, 'fric(\d+)', 'tokens', 'once');
    if ~isempty(tok), fric(i) = str2double(tok{1}); end
    tok = regexp(d(i).name, 'catchment(\d+)', 'tokens', 'once');
    if ~isempty(tok), catchment(i) = str2double(tok{1}); end
end

file = fullfile(mdir, cellstr(name));
T = table(file, name, region, fric, catchment, bytes);
T = sortrows(T, {'region','fric','name'});

end
