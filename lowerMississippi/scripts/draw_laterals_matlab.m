clc;
clear all;

%% Plot style
set(groot, ...
    'DefaultAxesFontName',        'Helvetica', ...
    'DefaultAxesFontSize',        18, ...
    'DefaultAxesXColor',          'k', ...
    'DefaultAxesYColor',          'k', ...
    'DefaultTextColor',           'k', ...
    'DefaultLegendFontSize',      18, ...
    'DefaultLegendTextColor',     'k', ...
    'DefaultLineLineWidth',       3,...
    'DefaultLineMarkerSize',      12);

%% Load data
tarbert = readmatrix('../segment0/geo/boundary_Q');
morganza = readmatrix('../segment0/geo/lateralDatas/qlat0');
bonnet = readmatrix('../segment0/geo/lateralDatas/qlat1');
bohemia = readmatrix('../segment0/geo/lateralDatas/qlat6');
headofpasses = readmatrix('../segment0/geo/boundary_h');

morganza(:,2) = -morganza(:,2);
bonnet(:,2)   = -bonnet(:,2);
bohemia(:,2)  = -bohemia(:,2);

t0 = datetime(2011,1,1);

%% Create figure and tiled layout
figure('Units','in','Position',[1 1 7 5]);

%% ===================== AX1 =====================
ax1 = gca();
hold(ax1,'on')

% Time for Tarbert
dt = seconds(tarbert(:,1));
time_datetime = t0 + dt;

% Left axis (Discharge)
yyaxis(ax1,'left')
h1 = plot(ax1,time_datetime, tarbert(:,2), ':', ...
    'Color','blue', ...
    'MarkerIndices',1:10:length(tarbert(:,2)),...
    'LineWidth',3.5);
ylabel(ax1,'Discharge (m^3/s)')

% Right axis (Flow depth) → ONLY on ax1
dt = seconds(headofpasses(:,1));
time_datetime = t0 + dt;

yyaxis(ax1,'right')
ylim(ax1,[20 25])
h5 = plot(ax1,time_datetime, headofpasses(:,2), '-', ...
    'Color','red', ...
    'MarkerIndices',1:16:length(headofpasses(:,2)),'LineWidth',2.5);
ylabel(ax1,'Flow depth (m)')

% Keep both y-axes black
ax1.YAxis(1).Color = 'black';
ax1.YAxis(2).Color = 'black';

% Separate legend for ax1
legend(ax1,[h1 h5], ...
    {'Tarbert Landing inflow', ...
     'Head of Passes'}, ...
    'Location','best', ...
    'FontSize',16);

box on

%% ===================== AX2 =====================
figure('Units','in','Position',[1 1 7 5]);
ax2 = gca();
hold(ax2,'on')
% title(ax1,'Upstream Inflow and Downstream Flow Depth')
% 
% title(ax2,'Lateral Outflows')
% force both plots to use a 10^4 scaling factor on the y-axis
ax1.YAxis(1).Exponent = 4;   % left y-axis of ax1
% ax1.YAxis(2).Exponent = 4;   % right y-axis of ax1
ax2.YAxis.Exponent     = 4;  % only left y-axis for ax2

% Morganza
dt = seconds(morganza(:,1));
time_datetime = t0 + dt;
h2 = plot(ax2,time_datetime, morganza(:,2), ':', ...
        'Color','blue', ...
    'MarkerIndices',1:10:length(morganza(:,2)));

% Bonnet Carré
dt = seconds(bonnet(:,1));
time_datetime = t0 + dt;
h3 = plot(ax2,time_datetime, bonnet(:,2), '-', ...
        'Color','red', ...
    'MarkerIndices',1:6:length(bonnet(:,2)));

% Bohemia
dt = seconds(bohemia(:,1));
time_datetime = t0 + dt;
h4 = plot(ax2,time_datetime, bohemia(:,2), '-.', ...
        'Color','black', ...
    'MarkerIndices',1:25:length(bohemia(:,2)));

ylabel(ax2,'Discharge (m^3/s)')

% Separate legend for ax2
legend(ax2,[h2 h3 h4], ...
    {'Morganza Spillway', ...
     'Bonnet Carre Spillway', ...
     'Bohemia Spillway'}, ...
    'Location','best', ...
    'FontSize',16);
box on

text(ax1, 0.02, 0.95, '(a)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',18)

text(ax2, 0.02, 0.95, '(b)', ...
    'Units','normalized', ...
    'FontWeight','bold', ...
    'FontSize',18)

exportgraphics(ax1,'ex3_inflow.pdf','ContentType','vector');
exportgraphics(ax2,'ex3_lateral.pdf','ContentType','vector');