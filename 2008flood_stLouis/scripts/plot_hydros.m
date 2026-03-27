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

saverton = readtable('saverton_inflow');
valleycity = readtable('valleycity_inflow');
hermann = readtable('Hermann_inflow');
chester = readtable('chester_stage');
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
fig = figure('Units','in','Position',[1 1 7 5]);
% --- Top: Baton Rouge (Upstream / Dynamic)
ax1 = gca(); hold(ax1,'on');
% extract numeric vectors from tables using {}
mi = 1:round(max(1, numel(usgs_grafton.Date)/70)):numel(usgs_grafton.Date);
plot(ax1, usgs_grafton.Date, usgs_grafton.Q_cms, ...
    'LineWidth',3.5, 'LineStyle',':', 'Color','b', ...
    'MarkerIndices', mi, 'DisplayName','USGS');
% % % title(ax1,'Mississippi River at Grafton');
ylabel(ax1,'Discharge (m^3/s)');
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
mih = [15 60 105 150 195 235 270 310 350];

plot(ax1, hec_grafton.Date, hec_grafton.Q_cms, ...
    'LineStyle','--', ...
    'Color','r', ...
    'Marker','s', ...
    'MarkerIndices', mih, ...
    'MarkerFaceColor','r', ...
    'DisplayName','Dynamic (HEC-RAS)');


%% Read meshless grafton

Meshless_grafton = make_Meshless_df('../solution-ydk/segment3/run', 0, start_date);
mi  = [35 80 125 170 215 255 290 330 365];
plot(ax1, Meshless_grafton{:,3}, Meshless_grafton{:,2}, ...
    'LineStyle','-', ...
    'Color','k', ...
    'Marker','o', ...
    'MarkerIndices', mi, ...
    'MarkerFaceColor','k', ...
    'DisplayName','Meshless');
ax1.YAxis(1).Exponent = 4;
box on

%% st charles
fig = figure('Units','in','Position',[1 1 7 5]);
ax2 = gca(); hold(ax2,'on');
ax2.YAxis(1).Exponent = 4;
box on
mi = 1:round(max(1, numel(usgs_stCharles.Date)/70)):numel(usgs_stCharles.Date);
plot(ax2, usgs_stCharles.Date, usgs_stCharles.Q_cms, ...
    'LineWidth',3.5, 'LineStyle',':', 'Color','b', ...
    'MarkerIndices', mi, 'DisplayName','USGS');
% % title(ax2,'Missouri River at St. Charles');
ylabel(ax2,'Discharge (m^3/s)'); 

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
mih = [15 60 105 150 190 230 270 315 355];
plot(ax2, hec_stCharles.Date, hec_stCharles.Q_cms, ...
    'LineStyle','--', ...
    'Color','r', ...
    'Marker','s', ...
    'MarkerIndices', mih, ...
    'MarkerFaceColor','r', ...
    'DisplayName','Dynamic (HEC-RAS)');

% opts = detectImportOptions(fullfile(hecPath, 'ex4_denser.txt'), 'NumHeaderLines', 10);
% hec = readtable(fullfile(hecPath, 'ex4_denser.txt'), opts);
% hec_stCharles = hec(hec{:,3} == 39912, :);
% dateh = datetime(hec_stCharles{:,4},'InputFormat','ddMMMyyyy');
% mih = 1:round(max(1, numel(dateh)/50)):numel(dateh);
% % Plot
% plot(ax2, dateh, hec_stCharles{:,6}, ...
%     'LineStyle','-', ...
%     'Color','r', ...
%     'Marker','s', ...
%     'MarkerIndices', mi, ...
%     'DisplayName','Dynamic');
% 
%% Read meshless st charles

Meshless_stCharles = make_Meshless_df('../solution-ydk/segment1/run', 24, start_date);
mi  = [35 80 125 170 210 250 290 335 370];
plot(ax2, Meshless_stCharles{:,3}, Meshless_stCharles{:,2}, ...
    'LineStyle','-', ...
    'Color','k', ...
    'Marker','o', ...
    'MarkerIndices', mi, ...
    'MarkerFaceColor','k', ...
    'DisplayName','Meshless');
ax2.YLim = [0 14000];

%% st louis
fig = figure('Units','in','Position',[1 1 7 5]);
ax3 = gca(); hold(ax3,'on');
ax3.YAxis(1).Exponent = 4;
ax3.YLim = [0 27000];
box on
mi = [1];
plot(ax3, usgs_stLouis.Date, usgs_stLouis.Q_cms, ...
    'LineWidth',3.5, 'LineStyle',':', 'Color','b', ...
    'MarkerIndices', mi, 'DisplayName','USGS');
% % title(ax3,'Mississippi River at St. Louis');
ylabel(ax3,'Discharge (m^3/s)'); 
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
mih = [15 60 105 150 190 235 275 315 355];
% Plot
plot(ax3, hec_stLouis.Date, hec_stLouis.Q_cms, ...
    'LineStyle','--', ...
    'Color','r', ...
    'Marker','s', ...
    'MarkerIndices', mih, ...
    'MarkerFaceColor','r', ...
    'DisplayName','Dynamic (HEC-RAS)');

% opts = detectImportOptions(fullfile(hecPath, 'ex4_denser.txt'), 'NumHeaderLines', 10);
% hec = readtable(fullfile(hecPath, 'ex4_denser.txt'), opts);
% hec_stLouis = hec(hec{:,3} == 111687, :);
% dateh = datetime(hec_stLouis{:,4},'InputFormat','ddMMMyyyy');
% mih = 1:round(max(1, numel(dateh)/50)):numel(dateh);
% % Plot
% plot(ax3, dateh, hec_stLouis{:,6}, ...
%     'LineStyle','-', ...
%     'Color','r', ...
%     'Marker','s', ...
%     'MarkerIndices', mi, ...
%     'DisplayName','Dynamic');

%% Read meshless st louis

Meshless_stLouis = make_Meshless_df('../solution-ydk/segment4/run', 4, start_date);
mi  = [35 80 125 170 210 255 295 335 370];
plot(ax3, Meshless_stLouis{:,3}, Meshless_stLouis{:,2}, ...
    'LineStyle','-', ...
    'Color','k', ...
    'Marker','o', ...
    'MarkerFaceColor','k', ...
    'MarkerIndices', mi, ...
    'DisplayName','Meshless');


%%
% Meshless_illinois_in = make_Meshless_df('../solution-ydk/segment2/run', 17, start_date);
% writetable(Meshless_illinois_in,'illinois_inflowbc.xlsx')
% Meshless_miss_in = make_Meshless_df('../solution-ydk/segment1/run', 21, start_date);
% writetable(Meshless_miss_in,'miss_inflowbc.xlsx')
%%
text(ax1, 0.02, 0.95, '(a)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',18)

% text(ax2, 0.02, 0.95, '(b)', ...
%     'Units','normalized', ...
%     'FontWeight','bold', ...
%     'FontSize',18)
% text(ax3, 0.02, 0.95, '(c)', ...
%     'Units','normalized', ...
%     'FontWeight','bold', ...
%     'FontSize',18)

lgd1 = legend(ax1,'Location','northeast');
lgd2 = legend(ax2,'Location','northeast');
lgd3 = legend(ax3,'Location','northeast');
lgd1.FontSize = 15;
lgd2.FontSize = 16;
lgd3.FontSize = 16;
exportgraphics(ax1,'ex4_1.pdf','ContentType','vector');
exportgraphics(ax2,'ex4_2.pdf','ContentType','vector');
exportgraphics(ax3,'ex4_3.pdf','ContentType','vector');

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
