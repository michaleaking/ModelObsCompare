function cfg = setup_moc(varargin)
%SETUP_MOC  Put ISSM and the ModelObsCompare toolkit on the MATLAB path.
%
%   Model vs Observation Glacier Velocity Comparison Module.
%
%   Run this once at the start of a session before calling any moc_* function.
%   For the ISSM *binary* distribution the path setup is simply bin/ + lib/
%   (there is no source tree / devpath to source).
%
%   Data locations are NOT hard-coded: each is resolved on first use by
%   MOC_PATH (session override -> environment variable -> shared config file ->
%   interactive prompt). Nothing is asked for here, so a session that only
%   touches observations never has to answer for the model folder. Run
%   MOC_CONFIGURE('-set') to fill everything in up front.
%
%   ISSM is optional: without it the observation, ROI and flowline tools all
%   work, and only the model-side functions (MOC_LOAD_MODEL and friends) fail.
%
%   Usage:
%      cfg = setup_moc();
%      cfg = setup_moc('quiet',true);    % no summary printout
%
%   Name/value options:
%      'quiet' : (logical) suppress the path summary (default false)
%
%   Output:
%      cfg : the configuration struct from CONFIG_MOC
%
%   See also CONFIG_MOC, MOC_PATH, MOC_CONFIGURE.

p = inputParser;
p.addParameter('quiet', false, @(v) islogical(v)||isnumeric(v));
p.parse(varargin{:});
quiet = logical(p.Results.quiet);

% --- this toolkit (recursively add all sub-folders) ---
here = fileparts(mfilename('fullpath'));
addpath(genpath(here));

% --- ISSM binary distribution (optional) ---
issmOk = false;
try
    issm_dir = moc_path('issm_dir', 'prompt', false);
    issm_bin = fullfile(issm_dir, 'bin');
    issm_lib = fullfile(issm_dir, 'lib');
    if isfolder(issm_bin)
        addpath(issm_bin);
        if isfolder(issm_lib)
            addpath(issm_lib);   % compiled mex files (InterpFromMeshToGrid_matlab, etc.)
        end
        issmOk = true;
    else
        warning('setup_moc:issmMissing', ...
            ['ISSM bin folder not found at %s — model-side functions will not ' ...
             'work. Fix it with moc_configure(''issm_dir'', ''/path/to/issm'').'], ...
            issm_bin);
    end
catch
    % issm_dir not configured: only worth mentioning in the summary
end

cfg = config_moc();

if ~quiet
    fprintf('\nModel vs Observation Glacier Velocity Comparison Module (moc) ready.\n');
    if issmOk
        fprintf('  ISSM             = %s\n', issm_dir);
    else
        fprintf('  ISSM             = (not configured — model-side functions disabled)\n');
    end
    S = moc_path_settings();
    unset = {};
    for k = 1:numel(S)
        n = S(k).name;
        if strcmp(n, 'issm_dir'), continue; end
        if isfield(cfg, n)
            fprintf('  %-16s = %s\n', n, cfg.(n));
        else
            unset{end+1} = n; %#ok<AGROW>
        end
    end
    if ~isempty(unset)
        fprintf('  not configured   : %s\n', strjoin(unset, ', '));
        fprintf('  -> run moc_configure(''-set'') to fill these in (asked once, then remembered)\n');
    end
end

end
