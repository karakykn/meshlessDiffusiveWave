clc;
clear all;

saverton = readtable('saverton_inflow');
valleycity = readtable('valleycity_inflow');
hermann = readtable('Hermann_inflow');
chester = readtable('chester_stage');


figure('Units','in','Position',[1 1 7 4]);

% ---- Left y-axis: Discharge ----
yyaxis left
x = saverton{:,2};
y = saverton{:,3};

h1 = plot(x, y, '-', ...
    'Color', 'black', ...
    'Marker', 'd', ...
    'MarkerIndices', 1:10:length(x));
hold on

x = valleycity{:,2};
y = valleycity{:,3};
h2 = plot(x, y, '-', ...
    'Color', 'blue', ...
    'Marker', 'o', ...
    'MarkerIndices', 1:10:length(x));
hold on

x = hermann{:,2};
y = hermann{:,3};
h3 = plot(x, y, '-', ...
    'Color', 'red', ...
    'Marker', 's', ...
    'MarkerIndices', 1:10:length(x));
hold on
ylabel('Discharge (cms)')

x = chester{:,2};
y = chester{:,3};
yyaxis right
h4 = plot(x, y, '--', ...
    'Color', 'black', ...
    'Marker', 'd', ...
    'MarkerIndices', 1:10:length(x), ...
    'LineWidth',.5);
hold on
ylabel('Flow depth (m)')

% ---- Fix axis + label colors ----
ax = gca;
ax.YAxis(1).Color = 'black';
ax.YAxis(2).Color = 'black';

% ---- Legend ----
legend([h1 h2 h3 h4 ], ...
       {'Saverton Inflow (Q)','Valley City Inflow (Q)','Hermann Inflow (Q)','Chester Flow Depth (h)'}, ...
       'Location','best');
