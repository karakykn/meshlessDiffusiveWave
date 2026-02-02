% ex2_compare_debug.m
clear; close all; clc;

%% Plot style
set(groot, 'DefaultAxesFontName','Helvetica', ...
    'DefaultAxesFontSize',12, ...
    'DefaultLineLineWidth',1.2);

%% Paths
caseName = '../';
hecras_dir = fullfile('..','data','ex2');
nwm_dir = fullfile('..','data','ex2_nwm');

fprintf('Working directory: %s\n', pwd);
fprintf('Checking expected data folders:\n');
fprintf('  HEC-RAS folder: %s -> %s\n', hecras_dir, booleanToString(isfolder(hecras_dir)));
fprintf('  NWM folder:     %s -> %s\n', nwm_dir, booleanToString(isfolder(nwm_dir)));

%% locate segment directories
segmentPattern = fullfile(caseName,'segment*');
segmentListing = dir(segmentPattern);
segmentNames = {segmentListing([segmentListing.isdir]).name};
fprintf('Found %d segment directories matching "%s"\n', numel(segmentNames), segmentPattern);

%% Load NWM (if present)
nwm_x = []; nwm_t = [];
x_file = fullfile(nwm_dir,'t_dsQ_CNX');
t_file = fullfile(nwm_dir,'t_dsQ_CNT');
if isfile(x_file)
    try
        nwm_x = readmatrix(x_file,'NumHeaderLines',1);
        fprintf('Loaded %s: size = [%d %d]\n', x_file, size(nwm_x,1), size(nwm_x,2));
    catch ME
        warning('Failed to read %s: %s', x_file, ME.message);
    end
else
    fprintf('NWM CNX file not found: %s\n', x_file);
end
if isfile(t_file)
    try
        nwm_t = readmatrix(t_file,'NumHeaderLines',1);
        fprintf('Loaded %s: size = [%d %d]\n', t_file, size(nwm_t,1), size(nwm_t,2));
    catch ME
        warning('Failed to read %s: %s', t_file, ME.message);
    end
else
    fprintf('NWM CNT file not found: %s\n', t_file);
end

%% Find how many channels to plot by scanning hecras files present
% assume files named us1.txt, ds1.txt, us2.txt, ds2.txt, ...
maxChannels = 6; % conservative upper bound; adjust if needed
availableChannels = [];
for ch = 1:maxChannels
    usFile = fullfile(hecras_dir, sprintf('us%d.txt', ch));
    dsFile = fullfile(hecras_dir, sprintf('ds%d.txt', ch));
    if isfile(usFile) || isfile(dsFile)
        availableChannels(end+1) = ch; %#ok<SAGROW>
    end
end
if isempty(availableChannels)
    warning('No HEC-RAS us*/ds* files detected in %s. Plot will show empty tiles.', hecras_dir);
    availableChannels = 1:3; % still create 3 tiles so user sees them
end
nChannels = numel(availableChannels);
fprintf('Will attempt to plot %d channel(s): %s\n', nChannels, mat2str(availableChannels));

%% Set up tiled figure
fig = figure('Units','in','Position',[1 1 9 6]);
tiledlayout(nChannels,1, 'TileSpacing','compact', 'Padding','compact');

% Keep track of legend labels already used to avoid duplicates
legendUsed = containers.Map();

for idx_i = 1:nChannels
    ch = availableChannels(idx_i);
    ax = nexttile; hold(ax,'on');

    % read HEC files for this channel (robustly)
    us_file = fullfile(hecras_dir, sprintf('us%d.txt', ch));
    ds_file = fullfile(hecras_dir, sprintf('ds%d.txt', ch));
    us_arr = []; ds_arr = [];
    if isfile(us_file)
        try
            us_arr = readmatrix(us_file);
            fprintf('Channel %d: read %s size [%d %d]\n', ch, us_file, size(us_arr,1), size(us_arr,2));
        catch ME
            warning('Channel %d: failed reading %s (%s)\n', ch, us_file, ME.message);
        end
    end
    if isfile(ds_file)
        try
            ds_arr = readmatrix(ds_file);
            fprintf('Channel %d: read %s size [%d %d]\n', ch, ds_file, size(ds_arr,1), size(ds_arr,2));
        catch ME
            warning('Channel %d: failed reading %s (%s)\n', ch, ds_file, ME.message);
        end
    end

    % Extract columns with safe fallback
    us_Q = []; ds_Q = []; timeHec = [];
    if ~isempty(us_arr)
        if size(us_arr,2) >= 5
            us_Q = us_arr(2:end,5);
        elseif size(us_arr,2) >= 4
            us_Q = us_arr(2:end,4);
        else
            us_Q = [];
        end
    end
    if ~isempty(ds_arr)
        if size(ds_arr,2) >= 5
            ds_Q = ds_arr(2:end,5);
        elseif size(ds_arr,2) >= 4
            ds_Q = ds_arr(2:end,4);
        else
            ds_Q = [];
        end
    end
    if ~isempty(ds_Q)
        timeHec = linspace(0,24,length(ds_Q)); % hours
    end

    % Plot HEC-RAS if non-empty
    if ~isempty(us_Q) && any(~isnan(us_Q))
        h1 = plot(ax, timeHec, us_Q, 'k--', 'LineWidth', 0.9);
        label = 'Upstream (HEC)';
        addLegendOnce(label, h1, legendUsed);
    else
        fprintf('Channel %d: upstream HEC empty or all NaN, skipping plot.\n', ch);
    end
    if ~isempty(ds_Q) && any(~isnan(ds_Q))
        h2 = plot(ax, timeHec, ds_Q, 'b-', 'LineWidth', 1.2, ...
            'Marker','s', 'MarkerIndices', 1:6:max(1,length(timeHec)));
        label = 'Downstream (HEC)';
        addLegendOnce(label, h2, legendUsed);
    else
        fprintf('Channel %d: downstream HEC empty or all NaN, skipping plot.\n', ch);
    end

    % Meshless results (iterate segment*/run/<time>/Q)
    meshPlotted = false;
    for s = 1:numel(segmentNames)
        segName = segmentNames{s};
        runPath = fullfile(caseName, segName, 'run');
        if ~isfolder(runPath); continue; end

        dirListing = dir(runPath);
        dirListing = dirListing([dirListing.isdir]);
        tvals = []; down_qs = [];
        for k = 1:numel(dirListing)
            namek = dirListing(k).name;
            if strcmp(namek,'.') || strcmp(namek,'..') || startsWith(namek,'.'); continue; end
            tnum = str2double(namek);
            if isnan(tnum); continue; end
            Qfile = fullfile(runPath, namek, 'Q');
            if isfile(Qfile)
                Qdata = load(Qfile);
                if ~isempty(Qdata)
                    tvals(end+1,1) = tnum; %#ok<SAGROW>
                    down_qs(end+1,1) = Qdata(end); %#ok<SAGROW>
                end
            end
        end
        if isempty(tvals); continue; end
        [tvals_s, sidx] = sort(tvals);
        down_qs = down_qs(sidx);
        % only plot if we have numeric values
        if any(~isnan(down_qs))
            hmesh = plot(ax, tvals_s/3600, down_qs, 'k-', 'LineWidth', 1.0, ...
                'Marker','o', 'MarkerIndices', 1:round(max(1,length(tvals_s)/20)):length(tvals_s));
            addLegendOnce('Downstream (Meshless)', hmesh, legendUsed);
            meshPlotted = true;
        end
    end
    if ~meshPlotted
        fprintf('Channel %d: no meshless results plotted (no segment run/Q files found)\n', ch);
    end

    % NWM plotting for channel (if present); try to pick correct columns like python
    if ~isempty(nwm_x)
        % python: nwm_x has time column at col 1, then per-channel columns
        % columns mapping is data dependent; attempt to find sensible columns
        col_base = 1; % time in col 1
        % try downstream columns similar to previous script: cols 2..n
        time_nwm = nwm_x(:,1)/3600; % seconds -> hours
        % attempt two columns per channel if present
        % choose a safe column index that exists
        if size(nwm_x,2) >= 3
            c1 = min(2, size(nwm_x,2));
            c2 = min(3, size(nwm_x,2));
            if any(~isnan(nwm_x(:,c1)))
                hnwm1 = plot(ax, time_nwm, nwm_x(:,c1), 'r-', 'Marker','^', ...
                    'MarkerIndices', 1:6:max(1,length(time_nwm)));
                addLegendOnce('NWM series 1', hnwm1, legendUsed);
            end
            if any(~isnan(nwm_x(:,c2)))
                hnwm2 = plot(ax, time_nwm, nwm_x(:,c2), 'Color',[1 0.5 0], 'Marker','v', ...
                    'MarkerIndices', 1:6:max(1,length(time_nwm)));
                addLegendOnce('NWM series 2', hnwm2, legendUsed);
            end
        end
    end

    % finalize tile
    xlim(ax, [0 24]);
    ylim(ax, [19 26]);
    title(ax, sprintf('Channel %d', ch));
    ylabel(ax, 'Discharge (cms)');
    grid(ax,'on');

    % If nothing was plotted in this axes, write text so it's visible
    hlines = findobj(ax,'Type','line');
    if isempty(hlines)
        text(0.5,0.5,'No data for this channel','Units','normalized', ...
            'HorizontalAlignment','center','FontSize',11,'Color',[0.5 0.5 0.5]);
    end

    drawnow; % force figure update so canvas isn't empty
end

% final xlabel
xlabel(tiledlayout(fig), 'Time (hr)');

% --- helper nested functions ---
function s = booleanToString(b)
    if b; s = 'YES'; else s = 'NO'; end
end

function addLegendOnce(label, handle, mapObj)
    if ~isKey(mapObj, label)
        try
            set(handle,'DisplayName',label);
        catch
            % ignore
        end
        mapObj(label) = true;
    else
        % avoid adding legend entry duplicate: clear DisplayName
        try
            set(handle,'DisplayName','');
        catch
        end
    end
end
