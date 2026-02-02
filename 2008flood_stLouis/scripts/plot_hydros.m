clc; clear all;

%% Plot style
set(groot, 'DefaultAxesFontName','Helvetica', ...
    'DefaultAxesFontSize',12);
start_date = datetime(2008, 1, 1, 0, 0, 0);
%% Read usgs grafton
usgs_dir = fullfile('..', 'data', 'hydros');
usgs_grafton = readtable(fullfile(usgs_dir, 'grafton_mississippi'), ...
                         'FileType', 'text', ...
                         'Delimiter', '\t', ...
                         'NumHeaderLines', 24);

% Convert the date column to datetime
% (column name '20d' must exist exactly as in the file)
usgs_grafton.Date = datetime(usgs_grafton.('Var3'));

% Convert discharge from cfs to cms
usgs_grafton.Q_cms = usgs_grafton.('Var4') / 35.31466621266132;

%% Read usgs st charles
usgs_dir = fullfile('..', 'data', 'hydros');
usgs_stCharles = readtable(fullfile(usgs_dir, 'stCharles_mississippi'), ...
                         'FileType', 'text', ...
                         'Delimiter', '\t', ...
                         'NumHeaderLines', 24);

% Convert the date column to datetime
% (column name '20d' must exist exactly as in the file)
usgs_stCharles.Date = datetime(usgs_stCharles.('Var3'));

% Convert discharge from cfs to cms
usgs_stCharles.Q_cms = usgs_stCharles.('Var4') / 35.31466621266132;
%% Read usgs st louis
usgs_dir = fullfile('..', 'data', 'hydros');
usgs_stLouis = readtable(fullfile(usgs_dir, 'stLouis_mississippi'), ...
                         'FileType', 'text', ...
                         'Delimiter', '\t', ...
                         'NumHeaderLines', 24);

% Convert the date column to datetime
% (column name '20d' must exist exactly as in the file)
usgs_stLouis.Date = datetime(usgs_stLouis.('x20d'));

% Convert discharge from cfs to cms
usgs_stLouis.Q_cms = usgs_stLouis.('x14n') / 35.31466621266132;

%% Figure usgs grafton
fig = figure('Units','in','Position',[1 1 7 12]);
tiled = tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

% --- Top: Baton Rouge (Upstream / Dynamic)
ax1 = nexttile; hold(ax1,'on');
% extract numeric vectors from tables using {}
mi = 1:round(max(1, numel(usgs_grafton.Date)/70)):numel(usgs_grafton.Date);
plot(ax1, usgs_grafton.Date, usgs_grafton.Q_cms, 'LineStyle','-', 'Color','b', ...
    'Marker','d', 'MarkerIndices', mi, 'DisplayName','USGS (Gage: 05587450)');
title(ax1,'Mississippi River at Grafton');
ylabel(ax1,'Discharge (cms)');
%% Read HEC grafton
hecPath = fullfile('..', 'data', 'hydros', 'ex4_hecresults');

% Read the CSV file, skipping the first 12 rows
hec_grafton = readtable(fullfile(hecPath, 'grafton.txt'), ...
                        'FileType', 'text', ...
                        'Delimiter', ',', ...
                        'NumHeaderLines', 12);

% Filter rows where the 2nd column equals 173029
hec_grafton = hec_grafton(hec_grafton.('Var2') == 173029, :);

% Convert the 3rd column to datetime (mixed formats)
rawDates = hec_grafton{:,3};

% If they are stored as strings or chars with quotes, remove quotes
rawDates = erase(rawDates, "'");

% Convert to datetime using explicit format
hec_grafton.Date = datetime(rawDates, ...
                            'InputFormat', 'ddMMMyyyy HHmm', ...
                            'Locale', 'en_US');

% Assign discharge column (4th column) to Q-cms
hec_grafton.Q_cms = hec_grafton.('Var4');
mi = 1:round(max(1, numel(hec_grafton.Date)/70)):numel(hec_grafton.Date);
% Plot
plot(ax1, hec_grafton.Date, hec_grafton.Q_cms, ...
    'LineStyle','-', ...
    'Color','r', ...
    'Marker','s', ...
    'MarkerIndices', mi, ...
    'DisplayName','Dynamic');

%% Read meshless grafton

Meshless_grafton = make_Meshless_df('../solution-ydk/segment3/run', 0, start_date);
mi = 1:round(max(1, numel(Meshless_grafton{:,3})/80)):numel(Meshless_grafton{:,3});
plot(ax1, Meshless_grafton{:,3}, Meshless_grafton{:,2}, ...
    'LineStyle','-', ...
    'Color','k', ...
    'Marker','o', ...
    'MarkerIndices', mi, ...
    'DisplayName','Meshless');

%% st charles
ax2 = nexttile; hold(ax2,'on');
mi = 1:round(max(1, numel(usgs_stCharles.Date)/70)):numel(usgs_stCharles.Date);
plot(ax2, usgs_stCharles.Date, usgs_stCharles.Q_cms, 'LineStyle','-', 'Color','b', ...
    'Marker','d', 'MarkerIndices', mi, 'DisplayName','USGS (Gage: 06935965)');
title(ax2,'Missouri River at St. Charles');
ylabel(ax2,'Discharge (cms)'); 

%% Read HEC st charles
hecPath = fullfile('..', 'data', 'hydros', 'ex4_hecresults');

% Read the CSV file, skipping the first 12 rows
hec_stCharles = readtable(fullfile(hecPath, 'stCharles.txt'), ...
                        'FileType', 'text', ...
                        'Delimiter', ',', ...
                        'NumHeaderLines', 12);

% Filter rows where the 2nd column equals 173029
hec_stCharles = hec_stCharles(hec_stCharles.('Var2') == 39912, :);

% Convert the 3rd column to datetime (mixed formats)
rawDates = hec_stCharles{:,3};

% If they are stored as strings or chars with quotes, remove quotes
rawDates = erase(rawDates, "'");

% Convert to datetime using explicit format
hec_stCharles.Date = datetime(rawDates, ...
                            'InputFormat', 'ddMMMyyyy HHmm', ...
                            'Locale', 'en_US');

% Assign discharge column (4th column) to Q-cms
hec_stCharles.Q_cms = hec_stCharles.('Var4');
mi = 1:round(max(1, numel(hec_stCharles.Date)/70)):numel(hec_stCharles.Date);
% Plot
plot(ax2, hec_stCharles.Date, hec_stCharles.Q_cms, ...
    'LineStyle','-', ...
    'Color','r', ...
    'Marker','s', ...
    'MarkerIndices', mi, ...
    'DisplayName','Dynamic');

%% Read meshless st charles

Meshless_stCharles = make_Meshless_df('../solution-ydk/segment1/run', 24, start_date);
mi = 1:round(max(1, numel(Meshless_stCharles{:,3})/80)):numel(Meshless_stCharles{:,3});
plot(ax2, Meshless_stCharles{:,3}, Meshless_stCharles{:,2}, ...
    'LineStyle','-', ...
    'Color','k', ...
    'Marker','o', ...
    'MarkerIndices', mi, ...
    'DisplayName','Meshless');

%% st louis
ax3 = nexttile; hold(ax3,'on');
mi = 1:round(max(1, numel(usgs_stLouis.Date)/70)):numel(usgs_stLouis.Date);
plot(ax3, usgs_stLouis.Date, usgs_stLouis.Q_cms, 'LineStyle','-', 'Color','b', ...
    'Marker','d', 'MarkerIndices', mi, 'DisplayName','USGS (Gage: 07010000)');
title(ax3,'Mississippi River at St. Louis');
ylabel(ax3,'Discharge (cms)'); 
%% Read HEC st louis
hecPath = fullfile('..', 'data', 'hydros', 'ex4_hecresults');

% Read the CSV file, skipping the first 12 rows
hec_stLouis = readtable(fullfile(hecPath, 'stLouis.txt'), ...
                        'FileType', 'text', ...
                        'Delimiter', ',', ...
                        'NumHeaderLines', 12);

% Filter rows where the 2nd column equals 173029
hec_stLouis = hec_stLouis(hec_stLouis.('Var2') == 111687, :);

% Convert the 3rd column to datetime (mixed formats)
rawDates = hec_stLouis{:,3};

% If they are stored as strings or chars with quotes, remove quotes
rawDates = erase(rawDates, "'");

% Convert to datetime using explicit format
hec_stLouis.Date = datetime(rawDates, ...
                            'InputFormat', 'ddMMMyyyy HHmm', ...
                            'Locale', 'en_US');

% Assign discharge column (4th column) to Q-cms
hec_stLouis.Q_cms = hec_stLouis.('Var4');
mi = 1:round(max(1, numel(hec_stLouis.Date)/70)):numel(hec_stLouis.Date);
% Plot
plot(ax3, hec_stLouis.Date, hec_stLouis.Q_cms, ...
    'LineStyle','-', ...
    'Color','r', ...
    'Marker','s', ...
    'MarkerIndices', mi, ...
    'DisplayName','Dynamic');
%% Read meshless st louis

Meshless_stLouis = make_Meshless_df('../solution-ydk/segment4/run', 4, start_date);
mi = 1:round(max(1, numel(Meshless_stLouis{:,3})/80)):numel(Meshless_stLouis{:,3});
plot(ax3, Meshless_stLouis{:,3}, Meshless_stLouis{:,2}, ...
    'LineStyle','-', ...
    'Color','k', ...
    'Marker','o', ...
    'MarkerIndices', mi, ...
    'DisplayName','Meshless');

% Add legends
legend(ax1,'Location','best');
legend(ax2,'Location','northwest');
legend(ax3,'Location','best');
% startX = datetime(2011,1,1);
% endX   = datetime(2012,1,1);
% 
% xlim(ax1, [startX endX]);
% xlim(ax2, [startX endX]);
%% Functions
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
