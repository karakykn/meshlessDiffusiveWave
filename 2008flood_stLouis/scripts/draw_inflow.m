clc;
clear all;
set(groot, ...
    'DefaultAxesFontName',        'Helvetica', ...
    'DefaultAxesFontSize',        16, ...
    'DefaultAxesXColor',          'k', ...
    'DefaultAxesYColor',          'k', ...
    'DefaultTextColor',           'k', ...
    'DefaultLegendFontSize',      12, ...
    'DefaultLegendTextColor',     'k', ...
    'DefaultLineLineWidth',       3.5,...
    'DefaultLineMarkerSize',      12);
saverton = readtable('saverton_inflow');
valleycity = readtable('valleycity_inflow');
hermann = readtable('Hermann_inflow');
chester = readtable('chester_stage');


figure('Units','in','Position',[1 1 7 5]);  % left y-axis of ax1
ax1 = gca();
hold(ax1,'on')
ax1.YAxis(1).Exponent = 4; 
% ax1.YAxis(2).Exponent = 4;   % right y-axis of ax1

% ---- Left y-axis: Discharge ----
yyaxis(ax1,'left')
x = saverton{:,2};
y = saverton{:,3};

h1 = semilogy(x, y, ':', ...
    'Color', 'blue', ...
    'MarkerIndices', 1:10:length(x));
hold on

x = valleycity{:,2};
y = valleycity{:,3};
h2 = semilogy(x, y, '-', ...
    'Color', 'red', ...
    'MarkerIndices', 1:10:length(x));
hold on

x = hermann{:,2};
y = hermann{:,3};
h3 = semilogy(x, y, '-.', ...
    'Color', 'k', ...
    'MarkerIndices', 1:10:length(x));
hold on
ylabel('Discharge (m^3/s)')

x = chester{:,2};
y = chester{:,3};
yyaxis(ax1,'right')
ylim(ax1,[-10 25])
h4 = semilogy(x, y, '--', ...
    'Color', 'm', ...
    'MarkerIndices', 1:10:length(x), ...
    'LineWidth',2.5);
hold on
ylabel('Flow depth (m)')

% ---- Fix axis + label colors ----
ax = gca;
ax.YAxis(1).Color = 'black';
ax.YAxis(2).Color = 'black';

% ---- Legend ----
lgd = legend([h1 h2 h3 h4 ], ...
       {'Saverton inflow','Valley City inflow','Hermann inflow','Chester'}, ...
       'Location','northwest');
lgd.FontSize = 16;
box on
% exportgraphics(ax,'ex4_inflows.pdf','ContentType','vector');