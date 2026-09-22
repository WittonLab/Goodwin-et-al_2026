%% ========================================================================
%  Extract imaging-window VR behaviour for GCaMP analysis
%  ========================================================================
%
%  Purpose
%  -------
%  This script crops Virmen behaviour variables to the behavioural epoch that
%  overlaps with a selected imaging session. It saves the cropped time,
%  position, velocity, reward and lick variables as ``*_tyv.mat`` files for
%  downstream GCaMP/Suite2p place-cell analysis.
%
%  Expected workspace inputs
%  -------------------------
%  This script assumes that the following variables have already been loaded
%  from the Virmen/behaviour loader:
%
%      vars    : Virmen metadata structure. Must contain vars.times.
%      t       : Behaviour time vector.
%      y       : VR position vector.
%      v       : Running velocity vector.
%      reward  : Reward vector, same length as t.
%      licks   : Lick vector, same length as t.
%
%  Output files
%  ------------
%  The script saves:
%
%      <mouse>_<date>_tyv.mat       full imaging-window behaviour
%      <mouse>_<date>_pre_tyv.mat   pre-switch behaviour segment
%      <mouse>_<date>_post_tyv.mat  post-switch behaviour segment
%
%  If two reward-position switches are found, a third post-switch epoch is
%  also saved:
%
%      <mouse>_<date>_post2_tyv.mat
%
%  Notes
%  -----
%  The original version of this script used hard-coded output filenames and
%  a fixed 3rd-to-4th laser timepoint rule. This cleaned version keeps that
%  behaviour as the default but makes it explicit and easier to change.

%% ========================================================================
%  1. User settings
%  ========================================================================

% Base output directory for cropped behaviour files.
% Edit this path for your machine/repository before running.
outputBaseDir = 'C:\Users\dg453\OneDrive - University of Exeter\Desktop\VR Analysis\Cohort 7\data';

% Which imaging epoch to extract.
% laserPairRows = [3 4] reproduces the previous behaviour: use the 3rd and
% 4th rows of vars.times.laser_times as imaging start and end.
laserPairRows = [3 4];

% Environment entries shorter than this are treated as artefacts/mistakes
% and ignored when identifying environment switches.
minEnvironmentDurationSec = 20;

% If true, a QC plot is displayed and you can manually choose the pre/post
% split point. This overrides automatic environment/reward switch splitting.
doManualSplit = false;

%% ========================================================================
%  2. Check required inputs
%  ========================================================================

requiredVars = {'vars','t','y','v','reward','licks'};
for iVar = 1:numel(requiredVars)
    if ~exist(requiredVars{iVar}, 'var')
        error('Required variable "%s" is missing from the workspace.', requiredVars{iVar});
    end
end

if ~isfield(vars, 'times') || ~isfield(vars.times, 'laser_times')
    error('vars.times.laser_times is required to identify the imaging window.');
end

%% ========================================================================
%  3. Infer mouse/date labels and create output folder
%  ========================================================================

% Try to infer mouse and date from vars.exper.exp_path. If this fails,
% manually define mouse and formatted_date in the fallback section below.
mouse = '';
formatted_date = '';

if isfield(vars, 'exper') && isfield(vars.exper, 'exp_path')
    exp_path = vars.exper.exp_path;

    if iscell(exp_path)
        exp_path = exp_path{1};
    end

    [~, tail] = fileparts(exp_path);
    mouseDateTokens = regexp(tail, '(DG\d+)_?(\d{6,8})', 'tokens');

    if ~isempty(mouseDateTokens)
        mouse = mouseDateTokens{1}{1};
        dateStr = mouseDateTokens{1}{2};

        if numel(dateStr) == 8
            formatted_date = datestr(datetime(dateStr, 'InputFormat', 'yyyyMMdd'), 'yymmdd');
        elseif numel(dateStr) == 6
            formatted_date = dateStr;
        end
    end
end

% Fallback labels if automatic parsing fails.
if isempty(mouse) || isempty(formatted_date)
    warning('Could not infer mouse/date from vars.exper.exp_path. Using fallback labels.');
    mouse = 'mouseID';
    formatted_date = datestr(datetime('today'), 'yymmdd');
end

outputDir = fullfile(outputBaseDir, formatted_date, mouse);
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

fprintf('Output folder: %s\n', outputDir);

%% ========================================================================
%  4. Extract behaviour during the imaging window
%  ========================================================================

laserTimes = vars.times.laser_times;

if isempty(laserTimes)
    error('No imaging laser times found in vars.times.laser_times.');
end

if size(laserTimes, 1) < max(laserPairRows)
    error('Requested laserPairRows [%d %d], but vars.times.laser_times only has %d rows.', ...
        laserPairRows(1), laserPairRows(2), size(laserTimes, 1));
end

tStart = laserTimes(laserPairRows(1), 2);
tEnd = laserTimes(laserPairRows(2), 2);

imagingIdx = find(t >= tStart & t <= tEnd);
if isempty(imagingIdx)
    error('No behaviour samples found between selected laser start/end times.');
end

st = imagingIdx(1);
ed = imagingIdx(end);
sessionDurationMin = (t(ed) - t(st)) / 60;

timg = t(imagingIdx);
yimg = y(imagingIdx);
vimg = v(imagingIdx);
rewardimg = reward(imagingIdx);
licksimg = licks(imagingIdx);

fprintf('Extracted imaging window %.2f to %.2f s: %.1f min.\n', ...
    tStart, tEnd, sessionDurationMin);

fullSessionFile = fullfile(outputDir, [mouse '_' formatted_date '_tyv.mat']);
save(fullSessionFile, 'timg', 'yimg', 'vimg', 'rewardimg', 'licksimg');
fprintf('Saved full imaging-window behaviour: %s\n', fullSessionFile);

%% ========================================================================
%  5. Identify environment switches within the session
%  ========================================================================

% Environment switches are useful for segmenting sessions with a change in
% the VR world/environment. Very short environments are removed because they
% often reflect accidental transitions.

envStartTimes = [];
envBounds = [];

if isfield(vars.times, 'currentWorld_times') && ~isempty(vars.times.currentWorld_times)
    envStartTimes = vars.times.currentWorld_times(:, 2);

    if numel(envStartTimes) > 1
        envDurations = diff(envStartTimes);
        keepEnv = [true; envDurations(:) > minEnvironmentDurationSec];
        envStartTimes = envStartTimes(keepEnv);
    end

    % Convert environment start/end times into full-session sample indices.
    for iEnv = 1:numel(envStartTimes)
        if iEnv == numel(envStartTimes)
            envIdx = find(t >= envStartTimes(iEnv));
        else
            envIdx = find(t >= envStartTimes(iEnv) & t <= envStartTimes(iEnv + 1));
        end

        if ~isempty(envIdx)
            envBounds(iEnv, :) = [envIdx(1), envIdx(end)]; %#ok<SAGROW>
        end
    end

    fprintf('%d valid environment epochs detected.\n', size(envBounds, 1));
else
    fprintf('No currentWorld_times field found; environment switch segmentation skipped.\n');
end

%% ========================================================================
%  6. Identify reward-position switches within the imaging window
%  ========================================================================

% Reward/object-position switches are used as the fallback split criterion
% when a clear environment transition is not available.

t_sw = [];

if isfield(vars.times, 'obj_change_times') && ~isempty(vars.times.obj_change_times)
    allRewardSwitchTimes = vars.times.obj_change_times(:, 2);
    t_sw = allRewardSwitchTimes(allRewardSwitchTimes >= min(timg) & allRewardSwitchTimes <= max(timg));
end

fprintf('%d reward-position switch(es) detected inside imaging window.\n', numel(t_sw));

%% ========================================================================
%  7. Split imaging-window behaviour into pre/post epochs
%  ========================================================================

% Priority order for splitting:
%   1. Manual split, if doManualSplit = true.
%   2. Environment switch, if there are at least two valid environments.
%   3. Reward-position switch, if one or two switches are detected.
%   4. Otherwise, split the imaging session in half for legacy compatibility.

splitMethod = '';

if doManualSplit
    figure;
    plot(timg, yimg, 'k');
    title('Click the desired pre/post split point');
    xlabel('Time (s)');
    ylabel('VR position');

    [xClick, ~] = ginput(1);
    [~, splitIdx] = min(abs(timg - xClick));
    splitMethod = 'manual';

elseif size(envBounds, 1) > 1
    % envBounds are in full t-space. Convert the first environment end to a
    % local index inside timg.
    splitIdx = envBounds(1, 2) - st + 1;
    splitIdx = max(1, min(splitIdx, numel(timg) - 1));
    splitMethod = 'environment';

elseif numel(t_sw) == 1
    [~, splitIdx] = min(abs(timg - t_sw(1)));
    splitMethod = 'single_reward_switch';

elseif numel(t_sw) >= 2
    [~, splitIdx] = min(abs(timg - t_sw(1)));
    [~, splitIdx2] = min(abs(timg - t_sw(2)));
    splitMethod = 'double_reward_switch';

else
    splitIdx = round(numel(timg) / 2);
    splitMethod = 'half_session';
    fprintf('No switch detected. Splitting imaging window in half.\n');
end

fprintf('Split method: %s\n', splitMethod);

% First pre/post split.
t_pre = timg(1:splitIdx);
y_pre = yimg(1:splitIdx);
v_pre = vimg(1:splitIdx);
reward_pre = rewardimg(1:splitIdx);
licks_pre = licksimg(1:splitIdx);

if strcmp(splitMethod, 'double_reward_switch')
    t_post = timg(splitIdx + 1:splitIdx2);
    y_post = yimg(splitIdx + 1:splitIdx2);
    v_post = vimg(splitIdx + 1:splitIdx2);
    reward_post = rewardimg(splitIdx + 1:splitIdx2);
    licks_post = licksimg(splitIdx + 1:splitIdx2);

    t_post2 = timg(splitIdx2 + 1:end);
    y_post2 = yimg(splitIdx2 + 1:end);
    v_post2 = vimg(splitIdx2 + 1:end);
    reward_post2 = rewardimg(splitIdx2 + 1:end);
    licks_post2 = licksimg(splitIdx2 + 1:end);
else
    t_post = timg(splitIdx + 1:end);
    y_post = yimg(splitIdx + 1:end);
    v_post = vimg(splitIdx + 1:end);
    reward_post = rewardimg(splitIdx + 1:end);
    licks_post = licksimg(splitIdx + 1:end);
end

%% ========================================================================
%  8. Save segmented behaviour files
%  ========================================================================

preFile = fullfile(outputDir, [mouse '_' formatted_date '_pre_tyv.mat']);
postFile = fullfile(outputDir, [mouse '_' formatted_date '_post_tyv.mat']);

save(preFile, 't_pre', 'y_pre', 'v_pre', 'reward_pre', 'licks_pre');
save(postFile, 't_post', 'y_post', 'v_post', 'reward_post', 'licks_post');

fprintf('Saved pre-switch behaviour:  %s\n', preFile);
fprintf('Saved post-switch behaviour: %s\n', postFile);

if exist('t_post2', 'var')
    post2File = fullfile(outputDir, [mouse '_' formatted_date '_post2_tyv.mat']);
    save(post2File, 't_post2', 'y_post2', 'v_post2', 'reward_post2', 'licks_post2');
    fprintf('Saved second post-switch behaviour: %s\n', post2File);
end

%% ========================================================================
%  9. Optional quick QC plot
%  ========================================================================

figure('Color', 'w');
plot(timg, yimg, 'k');
hold on;
xline(timg(splitIdx), '--r', 'Split 1');
if exist('splitIdx2', 'var')
    xline(timg(splitIdx2), '--b', 'Split 2');
end
xlabel('Time (s)');
ylabel('VR position');
title(sprintf('%s %s behaviour split: %s', mouse, formatted_date, splitMethod), ...
    'Interpreter', 'none');
box off;
