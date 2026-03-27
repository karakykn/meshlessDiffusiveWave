clear; close all; clc;


%% Case folder
caseName = '../';
%% -------------------------
% Global plotting parameters
%% -------------------------
set(groot, ...
    'DefaultAxesFontName',        'Helvetica', ...
    'DefaultAxesFontSize',        16, ...
    'DefaultAxesXColor',          'k', ...
    'DefaultAxesYColor',          'k', ...
    'DefaultTextColor',           'k', ...
    'DefaultLegendFontSize',      12, ...
    'DefaultLegendTextColor',     'k', ...
    'DefaultLineLineWidth',       3,...
    'DefaultLineMarkerSize',      12);

%% -------------------------
% Read HEC-RAS outputs
%% -------------------------
file = fullfile('..','data','ex1_hecras','ex1_hec_Qh.txt');

opts = detectImportOptions(file, 'NumHeaderLines', 10);
hec = readtable(file, opts);
hec.Q = hec{:,4};
hec.h = hec{:,6} - hec{:,5};
hec_5k = hec(hec{:,3} ~= 20000 & hec{:,3} ~= 0, :);
hec_5k.mins = linspace(0, 1440,145)';
hec_5k.hr = hec_5k.mins / 60 ;


mi_h= [1];
mi_c=[1];
mi_m=[1];

% date_h = datetime(hec_up{:,3},'InputFormat','ddMMMyyyy');
% timeHec = hec_up{:,4} ./ 60;
% % Plot
% plot(ax1, dateh, hec_grafton{:,6}, ...
%     'LineStyle','-', ...
%     'Color','r', ...
%     'Marker','s', ...
%     'MarkerIndices', mi, ...
%     'DisplayName','Dynamic');


%% Time axis for HEC-RAS
totalTime = 86400;                 % seconds
% timeHec = linspace(0,24,height(hec_up));  % hours
m = 8;
%% -------------------------
% Plot HEC-RAS
%% -------------------------
figure('Units','in','Position',[1 1 7 5]);
hold on;

plot(hec_5k.hr, hec_5k.h, ...
        'Color', 'r', ...
    'DisplayName','Dynamic (HEC-RAS)', ...
    'MarkerIndices', mi_h, ...
    'Marker','s', ...
    'LineStyle','--', ...
    'MarkerFaceColor','red');

%% -------------------------
% Read NWM files
%% -------------------------
% nwm_file_x = fullfile('..','data','ex1_nwm','t_dsQ_CNX');
% nwm_file_t = fullfile('..','data','ex1_nwm','t_dsQ_CNT');
% 
% nwm_x = readmatrix(nwm_file_x,'NumHeaderLines',1);
% nwm_t = readmatrix(nwm_file_t,'NumHeaderLines',1);
% 
% plot(nwm_x(:,1)/3600, nwm_x(:,2), ...
%     'Color', 'm', ...
%     'DisplayName','CNS (Beg \it{et. al.} 2023)', ...
%     'MarkerIndices', mi_c, ...
%     'Marker','^', ...
%     'LineStyle','-.', ...
%     'MarkerFaceColor','m');
% mi=[35 72 130];
% % plot(nwm_t(:,1)/3600, nwm_t(:,2), 'Color',[1 0.5 0], ...
%     % 'DisplayName','CNT - Downstream', 'Marker','v', 'MarkerFaceColor',[1 0.5 0],'MarkerIndices',mi_c2);

%% -------------------------
% Loop over Meshless segments
%% -------------------------
% segmentDirs = dir(fullfile(caseName,'segment*'));
start_date = datetime(2011, 1, 1, 0, 0, 0);
Meshless_A = make_Meshless_df('../segment0/run', 15, start_date);

plot(Meshless_A{:,1}/3600, Meshless_A{:,2}, ...
    'Color', 'k', ...
    'DisplayName','Meshless)', ...
    'MarkerIndices', mi_c, ...
    'Marker','o', ...
    'LineStyle','-', ...
    'MarkerFaceColor','k');

%% -------------------------
% Final formatting
%% -------------------------
xlim([0 24]);
% ylim([19.5 25.5]);
xlabel('Time (hr)');
ylabel('Flow depth (m)');
legend('Location','northeast');
box on
% exportgraphics(gcf,'ex1_Q.pdf','ContentType','vector');

function T = make_Meshless_df(run_dir, node_id, start_date)
    % Collect Q values
    [secs, dis] = collect_Q_values(run_dir, node_id);
    [secs, fd] = collect_h_values(run_dir, node_id);

    % Sort by seconds (safety, like Python code)
    [secs, order] = sort(secs);
    dis = dis(order);
    fd = fd(order);

    % Create table (similar to pandas DataFrame)
    T = table(secs, fd, dis, ...
        'VariableNames', {'seconds', 'flow_depth_m', 'discharge_cms'});

    % Create Date column
    T.Date = start_date + seconds(T.seconds);
end

function [times, values] = collect_Q_values(base_dir, node, verbose)
%COLLECT_Q_VALUES Collect values from 'Q' files in numeric subdirectories
%
% Returns sorted arrays (times, values).

    if nargin < 1 || isempty(base_dir)
        base_dir = 'run';
    end
    if nargin < 2
        node = -1;
    end
    if nargin < 3
        verbose = false;
    end

    % MATLAB is 1-indexed; Python code assumes 0-indexed
    if node >= 0
        node = node + 1;
    end

    % Gather numeric subdirectories
    d = dir(base_dir);
    numeric_dirs = [];

    for i = 1:length(d)
        if ~d(i).isdir || ismember(d(i).name, {'.', '..'})
            continue
        end

        t = str2double(d(i).name);
        if ~isnan(t)
            numeric_dirs = [numeric_dirs; t, i]; %#ok<AGROW>
        else
            if verbose
                fprintf('Skipping non-numeric directory: %s\n', d(i).name);
            end
        end
    end

    % Sort by numeric time
    numeric_dirs = sortrows(numeric_dirs, 1);

    times   = [];
    values  = [];
    skipped = {};

    for k = 1:size(numeric_dirs, 1)
        t = numeric_dirs(k, 1);
        dir_idx = numeric_dirs(k, 2);
        name = d(dir_idx).name;

        qpath = fullfile(base_dir, name, 'Q');

        if ~isfile(qpath)
            skipped(end+1,:) = {name, 'no Q file'}; %#ok<AGROW>
            if verbose
                fprintf('Skipping %s: Q file not found\n', name);
            end
            continue
        end

        try
            data = load(qpath);
        catch ME
            skipped(end+1,:) = {name, ['load error: ' ME.message]}; %#ok<AGROW>
            if verbose
                fprintf('Skipping %s: error loading Q: %s\n', name, ME.message);
            end
            continue
        end

        if isempty(data)
            skipped(end+1,:) = {name, 'empty Q'}; %#ok<AGROW>
            if verbose
                fprintf('Skipping %s: Q is empty\n', name);
            end
            continue
        end

        try
            if isvector(data)
                val = data(node);
            else
                val = data(node, end);
            end
        catch ME
            skipped(end+1,:) = {name, ['indexing error: ' ME.message]}; %#ok<AGROW>
            if verbose
                fprintf('Skipping %s: indexing error: %s\n', name, ME.message);
            end
            continue
        end

        times(end+1,1)  = t;      %#ok<AGROW>
        values(end+1,1) = val;    %#ok<AGROW>
    end

    if verbose && ~isempty(skipped)
        fprintf('Skipped %d entries. Example skips:\n', size(skipped,1));
        disp(skipped(1:min(5,end), :))
    end
end

function [times, values] = collect_h_values(base_dir, node, verbose)
%COLLECT_Q_VALUES Collect values from 'Q' files in numeric subdirectories
%
% Returns sorted arrays (times, values).

    if nargin < 1 || isempty(base_dir)
        base_dir = 'run';
    end
    if nargin < 2
        node = -1;
    end
    if nargin < 3
        verbose = false;
    end

    % MATLAB is 1-indexed; Python code assumes 0-indexed
    if node >= 0
        node = node + 1;
    end

    % Gather numeric subdirectories
    d = dir(base_dir);
    numeric_dirs = [];

    for i = 1:length(d)
        if ~d(i).isdir || ismember(d(i).name, {'.', '..'})
            continue
        end

        t = str2double(d(i).name);
        if ~isnan(t)
            numeric_dirs = [numeric_dirs; t, i]; %#ok<AGROW>
        else
            if verbose
                fprintf('Skipping non-numeric directory: %s\n', d(i).name);
            end
        end
    end

    % Sort by numeric time
    numeric_dirs = sortrows(numeric_dirs, 1);

    times   = [];
    values  = [];
    skipped = {};

    for k = 1:size(numeric_dirs, 1)
        t = numeric_dirs(k, 1);
        dir_idx = numeric_dirs(k, 2);
        name = d(dir_idx).name;

        qpath = fullfile(base_dir, name, 'h');

        if ~isfile(qpath)
            skipped(end+1,:) = {name, 'no h file'}; %#ok<AGROW>
            if verbose
                fprintf('Skipping %s: h file not found\n', name);
            end
            continue
        end

        try
            data = load(qpath);
        catch ME
            skipped(end+1,:) = {name, ['load error: ' ME.message]}; %#ok<AGROW>
            if verbose
                fprintf('Skipping %s: error loading h: %s\n', name, ME.message);
            end
            continue
        end

        if isempty(data)
            skipped(end+1,:) = {name, 'empty h'}; %#ok<AGROW>
            if verbose
                fprintf('Skipping %s: h is empty\n', name);
            end
            continue
        end

        try
            if isvector(data)
                val = data(node);
            else
                val = data(node, end);
            end
        catch ME
            skipped(end+1,:) = {name, ['indexing error: ' ME.message]}; %#ok<AGROW>
            if verbose
                fprintf('Skipping %s: indexing error: %s\n', name, ME.message);
            end
            continue
        end

        times(end+1,1)  = t;      %#ok<AGROW>
        values(end+1,1) = val;    %#ok<AGROW>
    end

    if verbose && ~isempty(skipped)
        fprintf('Skipped %d entries. Example skips:\n', size(skipped,1));
        disp(skipped(1:min(5,end), :))
    end
end

function result = trapezoidal(y,x)
    n = length(x);
    
    dx = zeros(n-1,1);
    dx(:) = x(2:end) - x(1:end-1);
    
    coef = zeros(n,1);
    coef(1:end-1) = dx(:);
    coef(end) = dx(end);
    coef(2:end-1) = coef(2:end-1) + dx(2:end);
    
    result = sum(y .* coef) / 2;
end

function rmseVal = rmse(exactY, exactX, approxY, approxX)
%RMSE   Root-mean-square error after interpolating approxY to exactX
%   rmseVal = rmse(exactY, exactX, approxY, approxX)
%
% Inputs:
%   exactY  - vector of "true" y-values
%   exactX  - vector of x-values corresponding to exactY
%   approxY - vector of approximate y-values (may be on different x-grid)
%   approxX - vector of x-values corresponding to approxY
%
% Output:
%   rmseVal - scalar RMSE value
%
% Notes:
%   - approxY is interpolated to the exactX positions using linear
%     interpolation. Points outside the range of approxX produce NaNs and
%     are excluded from the RMSE calculation (with a warning).
%   - Input vectors must be column or row vectors of matching lengths
%
% Example:
%   exactX = 0:0.1:1; exactY = sin(exactX);
%   approxX = 0:0.2:1; approxY = sin(approxX) + 0.05*randn(size(approxX));
%   r = rmse(exactY, exactX, approxY, approxX);

    % Basic validation
    if ~isvector(exactX) || ~isvector(exactY) || numel(exactX) ~= numel(exactY)
        error('exactX and exactY must be vectors of the same length.');
    end
    if ~isvector(approxX) || ~isvector(approxY) || numel(approxX) ~= numel(approxY)
        error('approxX and approxY must be vectors of the same length.');
    end

    % Ensure column vectors for consistent indexing
    exactX = exactX(:);
    exactY = exactY(:);
    approxX = approxX(:);
    approxY = approxY(:);

    % Interpolate approxY to the exactX grid
    approxY_interp = interp1(approxX, approxY, exactX, 'linear', NaN);

    % Handle points that could not be interpolated (NaNs)
    nanIdx = isnan(approxY_interp);
    if all(nanIdx)
        error('No overlap between approxX and exactX (all interpolated values are NaN).');
    elseif any(nanIdx)
        warning('%d/%d points fell outside approxX and were excluded from RMSE calculation.', ...
                sum(nanIdx), numel(nanIdx));
        approxY_interp = approxY_interp(~nanIdx);
        exactY = exactY(~nanIdx);
    end

    % Compute RMSE
    diffs = approxY_interp - exactY;
    rmseVal = sqrt(mean(diffs.^2));
end

function  MAPE = meanPercentageError(exactY, exactX, approxY, approxX)
%MEANPERCENTAGEERROR   Mean percentage error and mean absolute percentage error
%   [MPE, MAPE] = meanPercentageError(exactY, exactX, approxY, approxX)
%
% Outputs:
%   MPE  - mean percentage error in percent (can be positive or negative)
%   MAPE - mean absolute percentage error in percent (always >= 0)
%
% Same interpolation rules/notes as rmse.m:
%   - approxY is interpolated onto exactX
%   - points where exactY == 0 are excluded (and warned about) to avoid
%     division by zero
%
% Example:
%   [mpe, mape] = meanPercentageError(exactY, exactX, approxY, approxX);

    % Basic validation
    if ~isvector(exactX) || ~isvector(exactY) || numel(exactX) ~= numel(exactY)
        error('exactX and exactY must be vectors of the same length.');
    end
    if ~isvector(approxX) || ~isvector(approxY) || numel(approxX) ~= numel(approxY)
        error('approxX and approxY must be vectors of the same length.');
    end

    % Ensure column vectors
    exactX = exactX(:);
    exactY = exactY(:);
    approxX = approxX(:);
    approxY = approxY(:);

    % Interpolate approxY to exactX
    approxY_interp = interp1(approxX, approxY, exactX, 'linear', NaN);

    % Remove NaN-interpolated points
    nanIdx = isnan(approxY_interp);
    if all(nanIdx)
        error('No overlap between approxX and exactX (all interpolated values are NaN).');
    elseif any(nanIdx)
        warning('%d/%d points fell outside approxX and were excluded from calculation.', ...
                sum(nanIdx), numel(nanIdx));
        approxY_interp = approxY_interp(~nanIdx);
        exactY = exactY(~nanIdx);
    end

    % Exclude exactY == 0 to avoid divide-by-zero (report how many excluded)
    zeroIdx = exactY == 0;
    if all(zeroIdx)
        error('All exactY values are zero, cannot compute percentage errors.');
    elseif any(zeroIdx)
        warning('%d/%d points have exactY == 0 and were excluded from percentage calculation.', ...
                sum(zeroIdx), numel(zeroIdx));
        approxY_interp = approxY_interp(~zeroIdx);
        exactY = exactY(~zeroIdx);
    end

    % Percentage errors (in percent)
    pctErrors = (approxY_interp - exactY) ./ exactY * 100;

    MPE = mean(pctErrors);        % mean percentage error (signed)
    MAPE = mean(abs(pctErrors)); % mean absolute percentage error (non-negative)
end
