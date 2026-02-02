% ex2_compare_fixed2.m
clear; close all; clc;

%% Plot style
set(groot, 'DefaultAxesFontName','Helvetica', ...
    'DefaultAxesFontSize',12);

%% Paths
caseName = '../';
hecras_dir = fullfile('..','data','ex2');
nwm_dir = fullfile('..','data','ex2_nwm');

fprintf('HEC folder: %s\nNWM folder: %s\n', hecras_dir, nwm_dir);

%% Load NWM (assume first row header)
x_file = fullfile(nwm_dir,'t_dsQ_CNX');
t_file = fullfile(nwm_dir,'t_dsQ_CNT');

nwm_x = [];
nwm_t = [];
if isfile(x_file)
    nwm_x = readmatrix(x_file,'NumHeaderLines',1);
    fprintf('Loaded %s size [%d %d]\n', x_file, size(nwm_x,1), size(nwm_x,2));
else
    warning('Missing file: %s', x_file);
end
if isfile(t_file)
    nwm_t = readmatrix(t_file,'NumHeaderLines',1);
    fprintf('Loaded %s size [%d %d]\n', t_file, size(nwm_t,1), size(nwm_t,2));
else
    warning('Missing file: %s', t_file);
end

%% Determine available channels from hecras files
maxChannels = 6;
availableChannels = [];
for ch = 1:3   % you said 3 channels — keep it explicit
    if isfile(fullfile(hecras_dir, sprintf('us%d.txt',ch))) || ...
       isfile(fullfile(hecras_dir, sprintf('ds%d.txt',ch)))
        availableChannels(end+1) = ch; %#ok<SAGROW>
    end
end
if isempty(availableChannels)
    availableChannels = 1:3;
end
fprintf('Will plot channels: %s\n', mat2str(availableChannels));

%% Find segments
segmentListing = dir(fullfile(caseName,'segment*'));
segmentNames = {segmentListing([segmentListing.isdir]).name};
segmentNames = segmentNames(~ismember(segmentNames,{'.','..'}));
fprintf('Found segments: %s\n', strjoin(segmentNames, ', '));

%% Figure layout
nChannels = numel(availableChannels);
fig = figure('Units','in','Position',[1 1 7 12]);
tiled = tiledlayout(nChannels,1,'TileSpacing','compact','Padding','compact');

for ii = 1:nChannels
    ch = availableChannels(ii);
    ax = nexttile; hold(ax,'on');

    % --- Read HEC-RAS us/ch & ds/ch
    us_file = fullfile(hecras_dir, sprintf('us%d.txt', ch));
    ds_file = fullfile(hecras_dir, sprintf('ds%d.txt', ch));
    us_arr = []; ds_arr = [];
    if isfile(us_file); us_arr = readmatrix(us_file); end
    if isfile(ds_file); ds_arr = readmatrix(ds_file); end

    % robust extraction of Q column (4th or 5th)
    us_Q = []; ds_Q = []; timeHec = [];
    if ~isempty(us_arr)
        if size(us_arr,2) >= 5
            us_Q = us_arr(2:end,5);
        elseif size(us_arr,2) >= 4
            us_Q = us_arr(2:end,4);
        end
    end
    if ~isempty(ds_arr)
        if size(ds_arr,2) >= 5
            ds_Q = ds_arr(2:end,5);
        elseif size(ds_arr,2) >= 4
            ds_Q = ds_arr(2:end,4);
        end
    end
    if ~isempty(ds_Q)
        timeHec = linspace(0,24,length(ds_Q)); % hours
    end

    % Plot HEC upstream + downstream (try to match your single-case look)
    if ~isempty(us_Q) && any(~isnan(us_Q))
        if ii == 3
                   plot(ax, timeHec, us_Q, 'k--', 'Color', 'red', ...
             'Marker','s','MarkerIndices', 1:round(max(1,length(timeHec)/25)):length(timeHec), ...
            'DisplayName','Upstream (Dynamic)', 'LineWidth', .5);
        else
        plot(ax, timeHec, us_Q, 'k--', ...
             'MarkerIndices', 1:round(max(1,length(timeHec)/25)):length(timeHec), ...
            'DisplayName','Upstream', 'LineWidth', .5);
        end
    end
    if ~isempty(ds_Q) && any(~isnan(ds_Q))
        plot(ax, timeHec, ds_Q, 'r-', ...
            'Marker','s', 'MarkerIndices', 1:round(max(1,length(timeHec)/25)):length(timeHec), ...
            'DisplayName','Downstream (Dynamic)');
    end


    % --- NWM plots: explicit mapping for 7-column files
    % Expectation: col1=time, cols2-4=CNS for channels 1..3, cols5-7=CNT for channels 1..3
    if ~isempty(nwm_x)
        t_nwm = nwm_x(:,1)/3600; % seconds -> hours
        cns_col = 4 + ii;    % 2,3,4
        cnt_col = 4 + ii;     % 5,6,7
                if ii==3
                y_cns = nwm_x(:,4);
                step = max(1, round(length(t_nwm)/10));
                plot(ax, t_nwm, y_cns,'--', 'Color','m', ...
                    'Marker','^', 'MarkerIndices', 1:step:length(t_nwm), ...
                    'DisplayName', sprintf('Upstream (CNS)', cnt_col), 'LineWidth', .5);
                fprintf('Channel %d: plotting CNT col %d\n', ch, cnt_col);
        end
        if cns_col <= size(nwm_x,2)
            y_cns = nwm_x(:,cns_col);
            if any(~isnan(y_cns))
                step = max(1, round(length(t_nwm)/11));
                plot(ax, t_nwm, y_cns, 'm-', ...
                    'Marker','^', 'MarkerIndices', 1:step:length(t_nwm), ...
                    'DisplayName', sprintf('Downstream (CNS)', 'Interpreter', 'latex'));
                fprintf('Channel %d: plotting CNS col %d\n', ch, cns_col);
            else
                fprintf('Channel %d: CNS col %d all NaN\n', ch, cns_col);
            end
        else
            fprintf('Channel %d: CNS col %d out of range (ncols=%d)\n', ch, cns_col, size(nwm_x,2));
        end
        if ii==3
                y_cnt = nwm_x(:,4);
                step = max(1, round(length(t_nwm)/15));
                plot(ax, t_nwm, y_cnt,'--', 'Color',[1 0.5 0], ...
                    'Marker','v', 'MarkerIndices', 1:step:length(t_nwm), ...
                    'DisplayName', sprintf('Upstream (CNT)', cnt_col), 'LineWidth', .5);
                fprintf('Channel %d: plotting CNT col %d\n', ch, cnt_col);
        end

        if cnt_col <= size(nwm_x,2)
            y_cnt = nwm_x(:,cnt_col);
            
            if any(~isnan(y_cnt))
                step = max(1, round(length(t_nwm)/15));
                plot(ax, t_nwm, y_cnt, 'Color',[1 0.5 0], ...
                    'Marker','v', 'MarkerIndices', 1:step:length(t_nwm), ...
                    'DisplayName', sprintf('Downstream (CNT)', cnt_col));
                fprintf('Channel %d: plotting CNT col %d\n', ch, cnt_col);
            else
                fprintf('Channel %d: CNT col %d all NaN\n', ch, cnt_col);
            end
        else
            fprintf('Channel %d: CNT col %d out of range (ncols=%d)\n', ch, cnt_col, size(nwm_x,2));
        end

    end

    % --- Meshless results (plot downstream from each segment's run/<time>/Q)
    mesh_any = false;
    % for s = 1:length(segmentNames)
        runPath = fullfile(caseName, segmentNames{ii}, 'run');
        if ~isfolder(runPath); continue; end
        td = dir(runPath); td = td([td.isdir]);
        tvals = []; down_qs = []; up_qs = [];
        for k = 1:length(td)
            nm = td(k).name;
            if strcmp(nm,'.') || strcmp(nm,'..') || startsWith(nm,'.'); continue; end
            tnum = str2double(nm);
            if isnan(tnum); continue; end
            Qfile = fullfile(runPath, nm, 'Q');
            if isfile(Qfile)
                Qdata = load(Qfile);
                if ~isempty(Qdata)
                    tvals(end+1,1) = tnum; %#ok<SAGROW>
                    down_qs(end+1,1) = Qdata(end);
                    up_qs(end+1,1) = Qdata(1);%#ok<SAGROW>
                end
            end
        end

        if isempty(tvals); continue; end
        [tvals_s, sidx] = sort(tvals);
        down_qs = down_qs(sidx);
        up_qs = up_qs(sidx);
                if ii==3
            step = max(1, round(length(tvals_s)/20));
            plot(ax, tvals_s/3600, up_qs, 'k--', ...
                'Marker','o', 'MarkerIndices', 1:step:length(tvals_s), ...
                'DisplayName','Upstream (Meshless)', 'LineWidth', .5);
            mesh_any = true;
        end
        if any(~isnan(down_qs))
            step = max(1, round(length(tvals_s)/20));
            plot(ax, tvals_s/3600, down_qs, 'k-', ...
                'Marker','o', 'MarkerIndices', 1:step:length(tvals_s), ...
                'DisplayName','Downstream (Meshless)');
            mesh_any = true;
        end

    % end
    if ~mesh_any
        fprintf('Channel %d: no meshless runs found.\n', ch);
    end

    % finalize axes
    xlim(ax, [0 24]);

    % Auto y-limits from plotted lines (avoid fixed [19 26])
    lines = findobj(ax,'Type','line');
    if ~isempty(lines)
        yvals = [];
        for L = 1:length(lines)
            y = get(lines(L),'YData');
            y = y(~isnan(y) & isfinite(y));
            yvals = [yvals; y(:)]; %#ok<AGROW>
        end
        if ~isempty(yvals)
            ymin = min(yvals); ymax = max(yvals);
            pad = max(0.01*(ymax-ymin), 0.1);
            ylim(ax, [ymin-pad, ymax+pad]);
        else
            ylim(ax,'auto');
        end
    else
        % ylim(ax,'auto');
        text(0.5,0.5,'No data for this channel','Units','normalized','HorizontalAlignment','center');
    end

    title(ax, sprintf('Channel %d', ch));
    ylabel(ax, 'Discharge (cms)');
    legend(ax, 'Location','best');
    % if ii == 3
    %     ylim
    % grid(ax,'on');
    % if ii == 1
    %     legend(ax, 'Location','best', 'Interpreter','latex');
    % end
end

xlabel(tiled, 'Time (hr)');
legend(ax, 'Location','best');