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
file = fullfile('..','data','ex1_hecras','ex1_hec_n80.txt');

opts = detectImportOptions(file, 'NumHeaderLines', 10);
hec = readtable(file, opts);
hec_up = hec(hec{:,2} == 20000, :);
hec_down = hec(hec{:,2} == 0, :);


% mi_h= [10 27 40 55 80 110 140];
% mi_c=[15 32 47 55 67 90 125];
% mi_c2=[6 40 51 63 85 115];
% mi_m=[1 8 13 18 26 35 49];

mi_h= [17 31 61 112];
mi_c=[35 53 73 130];
% mi_c2=[49 85];
mi_m=[23 41 85];

% date_h = datetime(hec_up{:,3},'InputFormat','ddMMMyyyy');
% timeHec = hec_up{:,4} ./ 60;
% % Plot
% plot(ax1, dateh, hec_grafton{:,6}, ...
%     'LineStyle','-', ...
%     'Color','r', ...
%     'Marker','s', ...
%     'MarkerIndices', mi, ...
%     'DisplayName','Dynamic');

% ---- SAFETY: handle 4- or 5-column files ----
% Expected order (typical HEC-RAS):
% col 4 = stage, col 5 = discharge

%% Time axis for HEC-RAS
totalTime = 86400;                 % seconds
timeHec = linspace(0,24,height(hec_up));  % hours
m = 8;
%% -------------------------
% Plot HEC-RAS
%% -------------------------
figure('Units','in','Position',[1 1 7 5]);
hold on;

plot(timeHec, hec_up{:,3}, 'b', ...
    'DisplayName','Upstream inflow', 'LineStyle', ':', 'LineWidth',3.5);

plot(timeHec, hec_down{:,3}, ...
        'Color', 'r', ...
    'DisplayName','Dynamic (HEC-RAS)', ...
    'MarkerIndices', mi_h, ...
    'Marker','s', ...
    'LineStyle','--', ...
    'MarkerFaceColor','red');

%% -------------------------
% Read NWM files
%% -------------------------
nwm_file_x = fullfile('..','data','ex1_nwm','t_dsQ_CNX');
nwm_file_t = fullfile('..','data','ex1_nwm','t_dsQ_CNT');

nwm_x = readmatrix(nwm_file_x,'NumHeaderLines',1);
nwm_t = readmatrix(nwm_file_t,'NumHeaderLines',1);

plot(nwm_x(:,1)/3600, nwm_x(:,2), ...
    'Color', 'm', ...
    'DisplayName','CNS (Beg \it{et. al.} 2023)', ...
    'MarkerIndices', mi_c, ...
    'Marker','^', ...
    'LineStyle','-.', ...
    'MarkerFaceColor','m');
mi=[35 72 130];
% plot(nwm_t(:,1)/3600, nwm_t(:,2), 'Color',[1 0.5 0], ...
    % 'DisplayName','CNT - Downstream', 'Marker','v', 'MarkerFaceColor',[1 0.5 0],'MarkerIndices',mi_c2);

%% -------------------------
% Loop over Meshless segments
%% -------------------------
segmentDirs = dir(fullfile(caseName,'segment*'));

for i = 1:length(segmentDirs)
    if ~segmentDirs(i).isdir
        continue;
    end

    segmentPath = fullfile(caseName, segmentDirs(i).name, 'run');
    if ~isfolder(segmentPath)
        continue;
    end

    tvals = [];
    downQ = [];

    timeDirs = dir(segmentPath);
    for j = 1:length(timeDirs)
        if ~timeDirs(j).isdir || startsWith(timeDirs(j).name,'.')
            continue;
        end

        tnum = str2double(timeDirs(j).name);
        if isnan(tnum); continue; end

        Qfile = fullfile(segmentPath, timeDirs(j).name, 'Q');
        if isfile(Qfile)
            Qdata = load(Qfile);
            if ~isempty(Qdata)
                tvals(end+1,1) = tnum;
                downQ(end+1,1) = Qdata(end);
            end
        end
    end

    if isempty(tvals); continue; end

    [tvals, idx] = sort(tvals);
    downQ = downQ(idx);

    plot(tvals/3600, downQ,    'Color', 'k', ...
    'DisplayName','Meshless', ...
    'MarkerIndices', mi_m, ...
    'Marker','o', ...
    'MarkerFaceColor','k');
end

%% -------------------------
% Final formatting
%% -------------------------
xlim([0 24]);
ylim([19.5 25.5]);
xlabel('Time (hr)');
ylabel('Discharge (m^3/s)');
legend('Location','northeast');
box on
exportgraphics(gcf,'ex1_Q.pdf','ContentType','vector');