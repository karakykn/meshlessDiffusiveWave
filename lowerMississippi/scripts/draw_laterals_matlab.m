clc;
clear all;

tarbert = readmatrix('../segment0/geo/boundary_Q');
morganza = readmatrix('../segment0/geo/lateralDatas/qlat0');
bonnet = readmatrix('../segment0/geo/lateralDatas/qlat1');
bohemia = readmatrix('../segment0/geo/lateralDatas/qlat6');
headofpasses = readmatrix('../segment0/geo/boundary_h');
morganza(:,2) = -morganza(:,2);
bonnet(:,2) = -bonnet(:,2);
bohemia(:,2) = -bohemia(:,2);

t0 = datetime(2011,1,1);
dt = seconds(tarbert(:,1));
time_datetime = t0 + dt;

figure('Units','in','Position',[1 1 7 4]);

% ---- Left y-axis: Discharge ----
yyaxis left
h1 = plot(time_datetime, tarbert(:,2)/10, '-', ...
    'Color', 'black', ...
    'Marker', 'd', ...
    'MarkerIndices', 1:10:length(tarbert(:,2)));
ylabel('Discharge (cms)')
hold on

dt = seconds(morganza(:,1));
time_datetime = t0 + dt;
h2 = plot(time_datetime, morganza(:,2), '-', ...
    'Color', 'blue', ...
    'Marker', 'o', ...
    'MarkerIndices', 1:10:length(morganza(:,2)));

dt = seconds(bonnet(:,1));
time_datetime = t0 + dt;
h3 = plot(time_datetime, bonnet(:,2), '-', ...
    'Color', 'red', ...
    'Marker', 's', ...
    'MarkerIndices', 1:6:length(bonnet(:,2)));

dt = seconds(bohemia(:,1));
time_datetime = t0 + dt;
h4 = plot(time_datetime, bohemia(:,2), '-', ...
    'Color', 'm', ...
    'Marker', '^', ...
    'MarkerIndices', 1:25:length(bohemia(:,2)));

% ---- Right y-axis: Flow depth ----
dt = seconds(headofpasses(:,1));
time_datetime = t0 + dt;

yyaxis right
h5 = plot(time_datetime, headofpasses(:,2), '--', ...
    'LineWidth', .5, ...
    'Color', 'black', ...
    'Marker', 'd', ...
    'MarkerIndices', 1:16:length(headofpasses(:,2)));
ylabel('Flow depth (m)')

% ---- Fix axis + label colors ----
ax = gca;
ax.YAxis(1).Color = 'black';
ax.YAxis(2).Color = 'black';

% ---- Legend ----
legend([h1 h2 h3 h4 h5], ...
       {'Tarbert Landing Inflow (Q/10)','Morganza Spillway (-Q)','Bonnet Carré Spillway (-Q)','Bohemia Spillway (-Q)','Head of Passes (h)'}, ...
       'Location','best');
