clear; close all; clc;


%% Case folder
caseName = '../';
%% -------------------------
% Global plotting parameters
%% -------------------------
set(groot, ...
    'DefaultAxesFontName',        'Helvetica', ...
    'DefaultAxesFontSize',        18, ...
    'DefaultAxesXColor',          'k', ...
    'DefaultAxesYColor',          'k', ...
    'DefaultTextColor',           'k', ...
    'DefaultLegendFontSize',      20, ...
    'DefaultLegendTextColor',     'k', ...
    'DefaultLineLineWidth',       2.3,...
    'DefaultLineMarkerSize',      10);

%% -------------------------
% Read HEC-RAS outputs
%% -------------------------
file = fullfile('..','data','ex1_hecras','ex1hec.txt');

opts = detectImportOptions(file, 'NumHeaderLines', 10);
hec = readtable(file, opts);
hec_up = hec(hec{:,2} == 20000, :);
hec_down = hec(hec{:,2} == 0, :);


% mi_h= [10 27 40 55 80 110 140];
% mi_c=[15 32 47 55 67 90 125];
% mi_c2=[6 40 51 63 85 115];
% mi_m=[1 8 13 18 26 35 49];

% mi_h= [32 55 105];
% mi_c=[32 65 100];
% mi_c2=[49 85];
% mi_m=[8 16 26];

mi_h = 1:int32(height(hec)/10):height(hec);
mi_m = 4:int32(53/8):53;

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
figure('Units','in','Position',[1 2 7 5]);

%% ===================== AX1 =====================
ax = nexttile;
hold(ax,'on')

plot(ax, timeHec, hec_up{:,5}, 'k', ...
    'DisplayName','Upstream inflow', 'LineStyle', '-.', 'MarkerFaceColor','red');

h = plot(ax, timeHec, hec_down{:,5}, ...
    'DisplayName','Dynamic - Downstream', 'Marker','*', 'Color','red','MarkerIndices',mi_h, ...
    'LineStyle','--');

%% -------------------------
% Read NWM files
%% -------------------------
nwm_file_x = fullfile('..','data','ex1_nwm','t_dsQ_CNX');
nwm_file_t = fullfile('..','data','ex1_nwm','t_dsQ_CNT');

nwm_x = readmatrix(nwm_file_x,'NumHeaderLines',1);
nwm_t = readmatrix(nwm_file_t,'NumHeaderLines',1);

% plot(nwm_x(:,1)/3600, nwm_x(:,2), 'Color','m', ...
%     'DisplayName','CNS - Downstream', 'Marker','^', 'MarkerFaceColor','m','MarkerIndices',mi_c);

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

    h2 = plot(ax, tvals/3600, downQ, 'k', ...
        'DisplayName','Meshless Downstream', 'Marker','o','MarkerIndices',mi_m, ...
        'LineStyle','-');
end

%% -------------------------
% Final formatting
%% -------------------------
xlim([0 24]);
ylim([19.5 25.5]);
xlabel('Time (hr)');
ylabel('Discharge (m^3/s)');
legend('Location','northeast');
uistack(h, 'top')
box on