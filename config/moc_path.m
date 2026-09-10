function p = moc_path(name, varargin)
%MOC_PATH  Resolve one configured data location, asking for it if need be.
%
%   Nothing in this toolkit is hard-coded to one machine. A path is resolved
%   on demand, in this order:
%     1. an override set with MOC_CONFIGURE this session;
%     2. an environment variable (ISSM_DIR, MOC_MODEL_DIR, ...);
%     3. the shared config file (see MOC_CONFIG_FILE) — the Python side reads
%        and writes the same file, so configuring once covers both;
%     4. an interactive prompt, when MATLAB is not running with -batch; the
%        answer is saved to the config file, so each path is asked for once.
%
%   If none of those produce a path, the error names the setting and the three
%   ways to supply it.
%
%   Usage:
%      p = moc_path('model_dir')
%      p = moc_path('termini_shp', 'prompt',false)   % never ask; error instead
%      [p] = moc_path('roi_dir')                     % created if missing
%
%   Input:
%      name : (char) a setting from MOC_PATH_SETTINGS, e.g. 'obs_netcdf_dir'
%   Name/value options:
%      'prompt' : (logical) allow the interactive prompt (default true)
%
%   Output:
%      p : (char) the configured path
%
%   See also MOC_CONFIGURE, MOC_PATH_SETTINGS, CONFIG_MOC.

q = inputParser;
q.addParameter('prompt', true, @(v) islogical(v)||isnumeric(v));
q.parse(varargin{:});
allowPrompt = logical(q.Results.prompt);

name = char(name);
S = moc_path_settings();
k = find(strcmp({S.name}, name), 1);
if isempty(k)
    error('moc_path:unknown', 'Unknown path setting "%s"; known: %s.', ...
        name, strjoin({S.name}, ', '));
end
setting = S(k);

p = local_lookup(setting);

if isempty(p) && allowPrompt && local_is_interactive()
    p = local_ask(setting);
end

if isempty(p)
    error('moc_path:notConfigured', ...
        ['The path setting "%s" (%s) is not configured.\n' ...
         'Set it in any of these ways:\n' ...
         '  setenv(''%s'', ''/path/to/data'')\n' ...
         '  moc_configure(''%s'', ''/path/to/data'')\n' ...
         '  moc_configure(''-set'')     %% prompts, then remembers'], ...
        setting.name, setting.desc, setting.env, setting.name);
end

if strcmp(setting.kind, 'dir') && ~isfolder(p)
    if strcmp(setting.name, 'roi_dir')
        mkdir(p);
    end
end

end

% ===========================================================================
function p = local_lookup(setting)
%LOCAL_LOOKUP  Override -> environment -> config file -> built-in default.
p = '';

ov = moc_path_state('get');
if isfield(ov, setting.name) && ~isempty(ov.(setting.name))
    p = ov.(setting.name); return;
end

e = getenv(setting.env);
if ~isempty(e), p = moc_clean_path(e); return; end

[~, data] = moc_config_file();
if isfield(data, setting.name) && ~isempty(data.(setting.name))
    p = moc_clean_path(char(data.(setting.name))); return;
end

p = local_builtin(setting.name);
end

% ---------------------------------------------------------------------------
function p = local_builtin(name)
%LOCAL_BUILTIN  Defaults that ship with the repository.
p = '';
if strcmp(name, 'roi_dir')
    here = fileparts(fileparts(mfilename('fullpath')));   % config/ -> repo root
    p = fullfile(here, 'rois');
end
end

% ---------------------------------------------------------------------------
function tf = local_is_interactive()
%LOCAL_IS_INTERACTIVE  False under -batch, where INPUT cannot be answered.
try
    tf = ~batchStartupOptionUsed;
catch
    tf = true;    % older releases: assume a prompt can be answered
end
end

% ---------------------------------------------------------------------------
function p = local_ask(setting)
%LOCAL_ASK  Prompt for one path, validate it, and remember the answer.
p = '';
if strcmp(setting.kind, 'dir'), what = 'folder'; else, what = 'file'; end
fprintf('\nModelObsCompare needs the %s for "%s":\n    %s\n', ...
    what, setting.name, setting.desc);
fprintf('(press Enter to cancel; set %s to skip this next time)\n', setting.env);

for attempt = 1:3
    reply = input(sprintf('  path to %s: ', setting.name), 's');
    if isempty(strtrim(reply)), return; end
    cand = moc_clean_path(reply);
    ok = (strcmp(setting.kind,'dir') && isfolder(cand)) || ...
         (strcmp(setting.kind,'file') && isfile(cand));
    if ok
        moc_path_state('set', setting.name, cand);
        try
            f = moc_configure('-save');
            fprintf('  saved to %s\n', f);
        catch ME
            fprintf('  (could not save: %s)\n', ME.message);
        end
        p = cand;
        return;
    end
    fprintf('  no such %s: %s\n', what, cand);
end
end
