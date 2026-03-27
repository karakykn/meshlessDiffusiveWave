clc; clear all;

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


%% Read usgs baton
usgs_dir = fullfile('..','data','usgs_baton');
opts = detectImportOptions(usgs_dir,'Delimiter','\t');
usgs_baton = readtable(usgs_dir, opts);
usgs_baton(:,10) = usgs_baton(:,10) ./ 35.31466621266132;
usgs_baton(:,8) = usgs_baton(:,8) ./ 3.281;
dt_baton = datetime(usgs_baton{:,3}, 'InputFormat','yyyy-MM-dd');
%% Read usgs belle
usgs_dir = fullfile('..','data','usgs_belle');
opts = detectImportOptions(usgs_dir,'Delimiter','\t');
usgs_belle = readtable(usgs_dir, opts);
usgs_belle(:,6) = usgs_belle(:,6) ./ 35.31466621266132;
usgs_belle(:,4) = usgs_belle(:,4) ./ 3.281;
dt_belle = datetime(usgs_belle{:,3}, 'InputFormat','yyyy-MM-dd');
usgs_belle{:,6} = fillmissing(usgs_belle{:,6}, 'linear', 'EndValues','extrap');


%% Figure usgs baton
fig = figure('Units','in','Position',[1 1 7 5]);

% --- Top: Baton Rouge (Upstream / Dynamic)
ax1 = gca(); hold(ax1,'on');
% extract numeric vectors from tables using {}
x_baton = dt_baton;
y_baton = usgs_baton{:,10};         % numeric discharge (cms)
mi = 1:round(max(1, numel(x_baton)/5)):numel(x_baton);
plot(ax1, x_baton, y_baton, 'LineStyle',':', 'Color','b', ...
    'DisplayName','USGS','LineWidth',3.5);
% title(ax1,'Mississippi River at Baton Rouge');
ylabel(ax1,'Discharge (m^3/s)');
%% Read Hec-Ras
fname = fullfile('..','data', 'hec','lmr_hec.txt');
opts = detectImportOptions(fname, 'NumHeaderLines', 11);
hec = readtable(fname, opts);
hec = hec(hec{:,5} == 0, :);
hec_baton = hec(strcmp(hec{:,3}, hec{1,3}), :);
hec_belle = hec(strcmp(hec{:,3}, hec{400,3}), :);
dateh = datetime(hec_baton{:,4},'InputFormat','ddMMMyyyy');
mih = [75 149 270];
plot(ax1, dateh, hec_baton{:,6}, ...
    'LineStyle','--', ...
    'Color','r', ...
    'Marker','s', ...
    'MarkerIndices', mih, ...
    'MarkerFaceColor','red', ...
    'DisplayName','Dynamic (HEC-RAS)');

%% Read NWM baton
fname = fullfile('..','data', 'beg_data','baton_disc_beg');
CNX_baton = readtable(fname, 'Delimiter', '\t');

% Downsample
secs = CNX_baton{1:24:end, 1};   % use {} to get numeric array
q       = CNX_baton{1:24:end, 2};

% Define start date
t0 = datetime(2011, 1, 1);

% Convert seconds to datetime
date = t0 + seconds(secs);
mi = [235 290 408];
% Plot
plot(ax1, date, q, ...
    'LineStyle','-.', ...
    'Color','m', ...
    'Marker','^', ...
    'MarkerIndices', mi, ...
    'MarkerFaceColor','m', ...
    'DisplayName','CNS (Beg \it{et. al.} 2023)');

%% Read meshless baton

start_date = datetime(2011, 1, 1, 0, 0, 0);

% Meshless datasets
Meshless_baton = make_Meshless_df('../segment0/run', 28, start_date);
mi = [28 93 115 218];
plot(ax1, Meshless_baton{:,3}, Meshless_baton{:,2}, ...
    'LineStyle','-', ...
    'Color','k', ...
    'Marker','o', ...
    'MarkerFaceColor','k', ...
    'MarkerIndices', mi, ...
    'DisplayName','Meshless');
box on



%% Belle Chasse
fig = figure('Units','in','Position',[1 1 7 5]);

% --- Top: Baton Rouge (Upstream / Dynamic)
ax2 = gca(); hold(ax2,'on');
x_belle = dt_belle;
y_belle = usgs_belle{:,6};          % numeric discharge (cms)
mi2 = 1:round(max(1, numel(x_belle)/20)):numel(x_belle);
plot(ax2, x_belle, y_belle, 'LineStyle',':', 'Color','b', ...
    'MarkerIndices', mi, 'DisplayName','USGS',...
    'LineWidth',3.5);

mi = [88 130 180 350];
% title(ax2,'Mississippi River at Belle Chasse');
ylabel(ax2,'Discharge (m^3/s)');
plot(ax2, dateh, hec_belle{:,6}, ...
    'LineStyle','--', ...
    'Color','r', ...
    'Marker','s', ...
    'MarkerFaceColor','r', ...
    'MarkerIndices', mi, ...
    'DisplayName','Dynamic (HEC-RAS)');

%% Read NWM belle
fname = fullfile('..','data', 'beg_data','belle_disc_beg');
CNX_belle = readtable(fname, 'Delimiter', '\t');

% Downsample
secs = CNX_belle{1:24:end, 1};   % use {} to get numeric array
q       = CNX_belle{1:24:end, 2};

% Define start date
t0 = datetime(2011, 1, 1);

% Convert seconds to datetime
date = t0 + seconds(secs);
mi = [128 288 418 667];
% Plot
plot(ax2, date, q, ...
    'LineStyle','-.', ...
    'Color','m', ...
    'Marker','^', ...
    'MarkerFaceColor','m', ...
    'MarkerIndices', mi, ...
    'DisplayName','CNS (Beg \it{et. al.} 2023)');
% 
% 
%% Read meshless belle

Meshless_belle = make_Meshless_df('../segment0/run', 72, start_date);
mi = [47 108 170];
plot(ax2, Meshless_belle{:,3}, Meshless_belle{:,2}, ...
    'LineStyle','-', ...
    'Color','k', ...
    'Marker','o', ...
    'MarkerIndices', mi, ...
    'MarkerFaceColor','k', ...
    'DisplayName','Meshless');
% 
% 
%% Plot settings
% Improve date tick formatting (auto handles datetimes nicely)
ax1.XTickLabelRotation = 0;
ax2.XTickLabelRotation = 0;
% datetick(ax2,'x','dd-mmm-yyyy','keepticks'); % optional: format bottom ticks
% If you prefer automatic datetime tick formatting remove datetick line above.

% Add legends
legend(ax1, 'Interpreter', 'tex');
legend(ax2, 'Interpreter', 'tex');
startX = datetime(2011,1,1);
endX   = datetime(2012,1,1);

xlim(ax1, [startX endX]);
xlim(ax2, [startX endX]);
box on

text(ax1, 0.02, 0.95, '(a)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',18)

text(ax2, 0.02, 0.95, '(b)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',18)

exportgraphics(ax1,'ex3_1.pdf','ContentType','vector');
exportgraphics(ax2,'ex3_2.pdf','ContentType','vector');

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
