clear; close all; clc;

%% Plot style
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
%% Paths
caseName = '../';
hecras_dir = fullfile('..','data','ex2_hec_ab_n80.txt');
nwm_dir = fullfile('..','..','..','Junk/Model_Results/Y-Channel/CNX_Output/q.txt');

fprintf('HEC folder: %s\nNWM folder: %s\n', hecras_dir, nwm_dir);
start_date = datetime(2011, 1, 1, 0, 0, 0);

%% Load NWM (assume first row header)
x_file = fullfile(nwm_dir);
nwm_x = readmatrix(x_file);
nwm_x_c2 = nwm_x(nwm_x(:,2) == 2, :);
nwm_x_c3 = nwm_x(nwm_x(:,2) == 3, :);
nwm_A = nwm_x_c2(:, [1 12]);
nwm_B = nwm_x_c3(:, [1 18]);
nwm_A(:, 1) = nwm_A(:, 1) / 60;
nwm_B(:, 1) = nwm_B(:, 1) / 60;

%% HEC read
hec=readmatrix(hecras_dir, 'NumHeaderLines',12);
hec_A = hec(hec(:,2)==2,4);
hec_B = hec(hec(:,2)==3,4);

%% Meshless read
Meshless_A = make_Meshless_df('../segment1/run', 9,start_date);
Meshless_B = make_Meshless_df('../segment2/run', 15, start_date);
UP = readmatrix(fullfile('..','data','ex2/us1.txt'));

%% Plot
figure('Units','in','Position',[1 1 7 5]);
ax1 = gca(); hold(ax1,'on');

plot(ax1, nwm_A(:,1), UP(:,5), ...
    'Color', 'b', ...
    'DisplayName','Upstream inflow', ...
    'LineWidth',3.5, ...
    'LineStyle',':');
mi=[16 23 29 90];
plot(ax1, nwm_A(:,1), hec_A(:), ...
    'Color', 'r', ...
    'DisplayName','Dynamic (HEC-RAS)', ...
    'MarkerIndices', mi, ...
    'Marker','s', ...
    'LineStyle','--', ...
    'MarkerFaceColor','red');
mi=[29 44 120];
plot(ax1, nwm_A(:,1), nwm_A(:,2), ...
    'Color', 'm', ...
    'DisplayName','CNS (Beg \it{et. al.} 2023)', ...
    'MarkerIndices', mi, ...
    'Marker','^', ...
    'LineStyle','-.', ...
    'MarkerFaceColor','m');
mi=[18 30 58 152];
plot(ax1, Meshless_A{:,1}/3600, Meshless_A{:,2}, ...
    'Color', 'k', ...
    'DisplayName','Meshless', ...
    'MarkerIndices', mi, ...
    'Marker','o', ...
    'MarkerFaceColor','k');
box on


figure('Units','in','Position',[1 1 7 5]);
ax2 = gca(); hold(ax2,'on');
% plot(ax2, nwm_A(:,1), UP(:,5), ...
%     'Color', 'b', ...
%     'DisplayName','Upstream inflows', ...
%     'LineWidth',3.5, ...
%     'LineStyle',':');
mi=[21 29 56 114];
plot(ax2, nwm_B(:,1), hec_B(:), ...
    'Color', 'r', ...
    'DisplayName','Dynamic (HEC-RAS)', ...
    'MarkerIndices', mi, ...
    'Marker','s', ...
    'LineStyle','--', ...
    'MarkerFaceColor','red');
mi=[37 58 130];
plot(ax2, nwm_B(:,1), nwm_B(:,2), ...
    'Color', 'm', ...
    'DisplayName','CNS (Beg \it{et. al.} 2023)', ...
    'MarkerIndices', mi, ...
    'Marker','^', ...
    'LineStyle','-.', ...
    'MarkerFaceColor','m');
mi=[25 58 110];
plot(ax2, Meshless_B{:,1}/3600, Meshless_B{:,2}, ...
    'Color', 'k', ...
    'DisplayName','Meshless', ...
    'MarkerIndices', mi, ...
    'Marker','o', ...
    'MarkerFaceColor','k');

text(ax1, 0.02, 0.95, '(a)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',18)

text(ax2, 0.02, 0.95, '(b)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',18)

xlim(ax1,[0 24]);
xlim(ax2,[0 24]);
% % ylim([19.5 25.5]);
xlabel('Time (hr)');
ylabel('Discharge (m^3/s)');
legend(ax1,'Location','northeast', 'Interpreter','tex');
legend(ax2,'Location','northeast', 'Interpreter','tex');
box on
xlabel(ax1,'Time (hr)');
ylabel(ax1,'Discharge (m^3/s)');
xlabel(ax2,'Time (hr)');
ylabel(ax2,'Discharge (m^3/s)');
exportgraphics(ax1,'ex2_Q.pdf','ContentType','vector');
exportgraphics(ax2,'ex2_Q_2.pdf','ContentType','vector');

%% Error A
mass_in = 7200 * 242.5;
mass_out_d = trapezoidal(hec_A(:), nwm_A(:,1)*3600);
mbe_d = ((mass_in - mass_out_d) / mass_in) * 100
mass_out_c = trapezoidal(nwm_A(:,2), nwm_A(:,1)*3600);
mbe_c = ((mass_in - mass_out_c) / mass_in) * 100
mass_out_m = trapezoidal(Meshless_A{:,2}, Meshless_A{:,1});
mbe_m = ((mass_in - mass_out_m) / mass_in) * 100

rmse_cA = rmse(hec_A(:), nwm_A(:,1), nwm_A(:,2), nwm_A(:,1))
rmse_mA = rmse(hec_A(:), nwm_A(:,1)*3600, Meshless_A{:,2}, Meshless_A{:,1})

mpe_cA = meanPercentageError(hec_A(:), nwm_A(:,1), nwm_A(:,2), nwm_A(:,1))
mpe_mA = meanPercentageError(hec_A(:), nwm_A(:,1)*3600, Meshless_A{:,2}, Meshless_A{:,1})

%% Error B
mass_in = 7200 * 242.5 *2;
mass_out_d = trapezoidal(hec_B(:), nwm_A(:,1)*3600);
mbe_d = ((mass_in - mass_out_d) / mass_in) * 100
mass_out_c = trapezoidal(nwm_B(:,2), nwm_A(:,1)*3600);
mbe_c = ((mass_in - mass_out_c) / mass_in) * 100
mass_out_m = trapezoidal(Meshless_B{:,2}, Meshless_B{:,1});
mbe_m = ((mass_in - mass_out_m) / mass_in) * 100

rmse_cB = rmse(hec_B(:), nwm_A(:,1), nwm_B(:,2), nwm_A(:,1))
rmse_mB = rmse(hec_B(:), nwm_A(:,1)*3600, Meshless_B{:,2}, Meshless_B{:,1})

mpe_cB = meanPercentageError(hec_B(:), nwm_A(:,1), nwm_B(:,2), nwm_A(:,1))
mpe_mB = meanPercentageError(hec_B(:), nwm_A(:,1)*3600, Meshless_B{:,2}, Meshless_B{:,1})

%% Functs
function T = make_Meshless_df(run_dir, node_id, start_date)
    % Collect Q values
    [secs, dis] = collect_Q_values(run_dir, node_id);

    % Sort by seconds (safety, like Python code)
    [secs, order] = sort(secs);
    dis = dis(order);

    % Create table (similar to pandas DataFrame)
    T = table(secs, dis, ...
        'VariableNames', {'seconds', 'discharge_cms'});

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
