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

startDate = datetime(2008,6,1,0,0,0);
endDate   = datetime(2008,9,1,0,0,0);
figure('Units','in','Position',[1 1 8 5]);  % left y-axis of ax1
ax1 = gca();
hold(ax1,'on')

%%
quincy = readtable('../data/hydros/to_hec/flowcfs.xlsx');
n = height(quincy);

% Dates you want to skip
missingDates = datetime([2008 2008],[7 7],[21 22]);

% Build dates by walking forward and skipping missingDates
datesVec = datetime.empty(0,1);
cur = startDate;
while numel(datesVec) < n
    if cur > endDate
        error('Ran past endDate before generating %d dates. Increase endDate or reduce n.', n)
    end
    if ~ismember(cur, missingDates)
        datesVec(end+1,1) = cur; %#ok<SAGROW>
    end
    cur = cur + days(1);
end

% Assign to table
quincy.Date = datesVec;
% plot(ax1, quincy.Date, quincy{:,3})
%%
ilinoy = readtable('illinois_inflowbc.xlsx');
ilinoy = ilinoy(1:24:end-1,:);
ilinoy = ilinoy(~ismember(ilinoy{:,1}, missingDates), :);
plot(ax1, ilinoy{:,1}, ilinoy{:,2}+quincy{:,3},'DisplayName','Quincy (USACE) + Illinois (1d) inflows')
%%
misery = readtable('miss_inflowbc.xlsx');
% plot(ax1, misery{:,1}, misery{:,2})
%%
usgsGraf = readtable('../data/hydros/grafton_mississippi');
mask = usgsGraf{:,3} >= startDate & usgsGraf{:,3} <= endDate;
usgsGraf = usgsGraf(mask, :);
usgsGraf{:,4} = usgsGraf{:,4} / 35.31466621266132;
plot(ax1, usgsGraf{:,3}, usgsGraf{:,4}, 'DisplayName','Grafton (USGS)');
%%
hecGraf = readtable('../data/hydros/wse_grafton_hecco2.csv');
plot(ax1, hecGraf{:,3}, hecGraf{:,4}, 'DisplayName','Grafton (HEC-RAS 2D)');

box on
legend
exportgraphics(ax1,'misspLevees.pdf','ContentType','vector');