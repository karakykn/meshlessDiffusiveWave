% convert this matalb code to that:
% ex2_compare_fixed2.m
clear; close all; clc;

%% Plot style
set(groot, ...
    'DefaultAxesFontName',        'Helvetica', ...
    'DefaultAxesFontSize',        18, ...
    'DefaultAxesXColor',          'k', ...
    'DefaultAxesYColor',          'k', ...
    'DefaultTextColor',           'k', ...
    'DefaultLegendFontSize',      20, ...
    'DefaultLegendTextColor',     'k', ...
    'DefaultLineLineWidth',       2.5,...
    'DefaultLineMarkerSize',      10);
%% Paths
caseName = '../';
hecras_dir = fullfile('..','data','ex2_hec.txt');
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
% maxChannels = 6;
% availableChannels = [];
% for ch = 1:3   % you said 3 channels — keep it explicit
%     if isfile(fullfile(hecras_dir, sprintf('us%d.txt',ch))) || ...
%             isfile(fullfile(hecras_dir, sprintf('ds%d.txt',ch)))
%         availableChannels(end+1) = ch; %#ok<SAGROW>
%     end
% end
% if isempty(availableChannels)
%     availableChannels = 1:3;
% end
% fprintf('Will plot channels: %s\n', mat2str(availableChannels));

%% Find segments
segmentListing = dir(fullfile(caseName,'segment*'));
segmentNames = {segmentListing([segmentListing.isdir]).name};
segmentNames = segmentNames(~ismember(segmentNames,{'.','..'}));
fprintf('Found segments: %s\n', strjoin(segmentNames, ', '));

%% Figure layout
% -----------------------
% Replace the original "Figure layout" block with the code below
% -----------------------
h = 8;
c = 8;
m = 9;
nChannels = 3;
letters = {'(a)','(b)','(c)'}; % labels for LaTeX composition; extend if >3 channels
segosU = [25000 30000 20000];
segosD = [20000 20000 0];

opts = detectImportOptions(hecras_dir, 'NumHeaderLines', 10);
hec = readtable(hecras_dir, opts);

% mi_h= [[8 18 20 23 30 55 95 140];[9 18 22 29 50 80 110 140];[17 17 17 17 32 50 95 140]];
% mi_c=[[18 27 37 72 110];[25 35 52 90 130];[32 48 64 82 82]];
% mi_c2=[[22 31 40 85 120];[30 42 62 85 125];[43 51 73 122 122]];
% mi_m=[[1 27 39 46 65 105 150 220];[1 25 41 48 69 120 160 228];[14 14 14 14 14 72 150 220]];
% 
% mi_hu= [17 22 30 50 101];
% mi_cu=[22 30 47 80 110];
% mi_c2u=[24 36 52 97];
% mi_mu=[10 25 44 73 120 200];
% % mi_m=[];

mi_h= [[19 39 120];[27 53 120];[33 59 105]];
mi_c=[[20 34 70];[30 54 82];[44 70 110]];
mi_c2=[[27 70];[42 82];[1 1]];
mi_m=[[32 90 160];[42 77 160];[71 71 155]];

mi_hu= [17 22 30 50 101];
mi_cu=[22 30 47 90];
mi_c2u=[24 36 52 97];
mi_mu=[25 44 73 120 200];
% mi_m=[];

% Loop over channels and create ONE figure per channel
for ii = 1:nChannels
    fig = figure('Units','in','Position',[1 1 7 5]); % wide-ish single-panel figure
    ax = gca(); hold(ax,'on');
    heco = hec(hec{:,1} == ii, :);
    hec_up = heco(heco{:,2} == segosU(ii), :);
    hec_down = heco(heco{:,2} == segosD(ii), :);
    us_Q = hec_up{:,5};
    ds_Q = hec_down{:,5};
    timeHec = linspace(0,24, height(hec_up));
    ch = ii;

    % Plot HEC upstream + downstream

    if ii == 3
        plot(ax, timeHec, us_Q, 'k--', 'Color', 'red', ...
            'Marker','s','MarkerIndices', mi_hu, ...
            'DisplayName','Dynamic - Upstream');
    else
        plot(ax, timeHec, us_Q, 'k--', ...
            'MarkerIndices', 1:round(max(1,length(timeHec)/h)):length(timeHec), ...
            'DisplayName','Upstream inflow');
    end

    plot(ax, timeHec, ds_Q, 'r-', ...
        'Marker','s', 'MarkerIndices', mi_h(ii,:), ...
        'DisplayName','Dynamic - Downstream', 'MarkerFaceColor','red');


    % --- NWM plots
    if ~isempty(nwm_x)
        t_nwm = nwm_x(:,1)/3600; % seconds -> hours
        cns_col = 4 + ii;    % expected columns mapping (tweak if needed)
        cnt_col = 4 + ii;
        if ii==3
            y_cns = nwm_x(:,4);
            step = round(max(1,length(t_nwm)/25));
            plot(ax, t_nwm, y_cns,'--', 'Color','m', ...
                'Marker','^', 'MarkerIndices', mi_cu, ...
                'DisplayName', sprintf('CNS - Upstream'));
            fprintf('Channel %d: plotting CNT col %d\n', ch, cnt_col);
        end
        if cns_col <= size(nwm_x,2)
            y_cns = nwm_x(:,cns_col);
            if any(~isnan(y_cns))
                step = max(1, round(length(t_nwm)/21));
                plot(ax, t_nwm, y_cns, 'Color', 'm',...
                    'Marker','^', 'MarkerIndices', mi_c(ii,:), ...
                    'DisplayName', 'CNS - Downstream',  'MarkerFaceColor','m');
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
            % plot(ax, t_nwm, y_cnt,'--', 'Color',[1 0.5 0], ...
            %     'Marker','v', 'MarkerIndices', mi_c2u, ...
            %     'DisplayName','CNT - Upstream');
            fprintf('Channel %d: plotting CNT col %d\n', ch, cnt_col);
        end
        if cnt_col <= size(nwm_x,2)
            y_cnt = nwm_x(:,cnt_col);
            if any(~isnan(y_cnt))
                step = max(1, round(length(t_nwm)/15));
                % plot(ax, t_nwm, y_cnt, 'Color',[1 0.5 0], ...
                %     'Marker','v', 'MarkerIndices', mi_c2(ii,:), ...
                %     'DisplayName','CNT - Downstream',  'MarkerFaceColor',[1 0.5 0]);
                fprintf('Channel %d: plotting CNT col %d\n', ch, cnt_col);
            else
                fprintf('Channel %d: CNT col %d all NaN\n', ch, cnt_col);
            end
        else
            fprintf('Channel %d: CNT col %d out of range (ncols=%d)\n', ch, cnt_col, size(nwm_x,2));
        end
    end

    % --- Meshless results (downstream from each segment's run/<time>/Q)
    mesh_any = false;
    runPath = fullfile(caseName, segmentNames{ii}, 'run');
    if isfolder(runPath)
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
                    down_qs(end+1,1) = Qdata(end); %#ok<SAGROW>
                    up_qs(end+1,1) = Qdata(1); %#ok<SAGROW>
                end
            end
        end
        if ~isempty(tvals)
            [tvals_s, sidx] = sort(tvals);
            down_qs = down_qs(sidx);
            up_qs = up_qs(sidx);
            if ii==3
                step = max(1, round(length(tvals_s)/20));
                plot(ax, tvals_s/3600, up_qs, 'k--', ...
                    'Marker','o', 'MarkerIndices', mi_mu, ...
                    'DisplayName','Meshless - Upstream');
                mesh_any = true;
            end
            if any(~isnan(down_qs))
                step = max(1, round(length(tvals_s)/20));
                plot(ax, tvals_s/3600, down_qs, 'k-', ...
                    'Marker','o', 'MarkerIndices', mi_m(ii,:), ...
                    'DisplayName','Meshless - Downstream', 'MarkerFaceColor','k');
                mesh_any = true;
            end
        end
    else
        fprintf('Channel %d: run folder not found (%s)\n', ch, runPath);
    end

    if ~mesh_any
        fprintf('Channel %d: no meshless runs found.\n', ch);
    end

    % finalize axes
    xlim(ax, [0 24]);

    % Auto y-limits from plotted lines
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
        text(0.5,0.5,'No data for this channel','Units','normalized','HorizontalAlignment','center');
    end

    ylabel(ax, 'Discharge (m^3/s)');
    xlabel(ax, 'Time (hr)');
    legend(ax, 'Location','northeast');
    if ii == 3
        legend(ax, 'Location','northeast', 'FontSize', 13.5);
    end

    % Add small panel letter for LaTeX composition (top-left)
    labelStr = '(?)';
    if ii <= numel(letters)
        labelStr = letters{ii};
    else
        labelStr = sprintf('(%d)', ii);
    end
    % text(ax, 0.02, 0.95, labelStr, 'Units','normalized', ...
    %     'FontSize', 18, 'FontWeight','bold', 'HorizontalAlignment','left');
    box on;
    drawnow;

    % Save figure: vector PDF + PNG backup
    outname_base = sprintf('channel_%d', ch);
    try
        exportgraphics(fig, [outname_base, '.pdf'], 'ContentType','vector');
        exportgraphics(fig, [outname_base, '.png'], 'Resolution',300);
        fprintf('Saved %s.pdf and %s.png\n', outname_base, outname_base);
    catch ME
        % fallback for older MATLAB
        warning('exportgraphics failed: %s. Trying print fallback.', ME.message);
        print(fig, '-dpdf', [outname_base, '.pdf']);
        print(fig, '-dpng', '-r300', [outname_base, '.png']);
    end

    % optionally close to save memory:
    % close(fig);
end