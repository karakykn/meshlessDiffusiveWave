clear; close all; clc;


%% Case folder
caseName = '../';
set(groot, 'DefaultAxesFontName','Helvetica', ...
    'DefaultAxesFontSize',12);
%% -------------------------
% Read HEC-RAS outputs
%% -------------------------
us1_file = fullfile('..','data','ex1_hecras','us1.txt');
ds1_file = fullfile('..','data','ex1_hecras','ds1.txt');

us_arr = readmatrix(us1_file);
ds_arr = readmatrix(ds1_file);

% ---- SAFETY: handle 4- or 5-column files ----
% Expected order (typical HEC-RAS):
% col 4 = stage, col 5 = discharge
if size(us_arr,2) >= 5
    us_Q = us_arr(2:end,5);
    us_h = us_arr(2:end,4) - 2;
else
    us_Q = us_arr(2:end,4);
    us_h = us_arr(2:end,3) - 2;
end

if size(ds_arr,2) >= 5
    ds_Q = ds_arr(2:end,5);
    ds_h = ds_arr(2:end,4);
else
    ds_Q = ds_arr(2:end,4);
    ds_h = ds_arr(2:end,3);
end

%% Time axis for HEC-RAS
totalTime = 86400;                 % seconds
timeHec = linspace(0,24,length(ds_Q));  % hours

%% -------------------------
% Plot HEC-RAS
%% -------------------------
figure('Units','in','Position',[1 1 7 4]);
hold on;

plot(timeHec, us_Q, 'k', ...
    'DisplayName','Upstream', 'LineStyle', '--', 'LineWidth', .5);

plot(timeHec, ds_Q, 'b', ...
    'DisplayName','Downstream (Dynamic)', 'Marker','s', 'Color','red','MarkerIndices',1:6:length(timeHec), 'LineWidth', 1.2);

%% -------------------------
% Read NWM files
%% -------------------------
nwm_file_x = fullfile('..','data','ex1_nwm','t_dsQ_CNX');
nwm_file_t = fullfile('..','data','ex1_nwm','t_dsQ_CNT');

nwm_x = readmatrix(nwm_file_x,'NumHeaderLines',1);
nwm_t = readmatrix(nwm_file_t,'NumHeaderLines',1);

plot(nwm_x(:,1)/3600, nwm_x(:,2), 'Color','m', ...
    'DisplayName','Downstream (CNS - Beg et al. 2023)', 'Marker','^','MarkerIndices',1:6:length(timeHec), 'LineWidth', 1.2);

plot(nwm_t(:,1)/3600, nwm_t(:,2), 'Color',[1 0.5 0], ...
    'DisplayName','Downstream (CNT - Beg et al. 2023)', 'Marker','v','MarkerIndices',1:4:length(timeHec), 'LineWidth', 1.2);

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

    plot(tvals/3600, downQ, 'k', ...
        'DisplayName','Downstream (Meshless)', 'Marker','o','MarkerIndices',1:3:length(timeHec), 'LineWidth', 1.2);
end

%% -------------------------
% Final formatting
%% -------------------------
xlim([0 24]);
ylim([19 26]);
xlabel('Time (hr)');
ylabel('Discharge (cms)');
legend('Location','best');
