function h = moc_home_dir()
%MOC_HOME_DIR  The current user's home folder, without needing Java.
%
%   Uses the platform environment variable, so it works under -nodisplay and
%   -nojvm where java.lang.System would warn or be unavailable.
%
%   Usage:
%      h = moc_home_dir()
%
%   Output:
%      h : (char) home folder ('' if the environment does not say)

if ispc
    h = getenv('USERPROFILE');
    if isempty(h)
        h = fullfile(getenv('HOMEDRIVE'), getenv('HOMEPATH'));
    end
else
    h = getenv('HOME');
end
h = char(h);
end
