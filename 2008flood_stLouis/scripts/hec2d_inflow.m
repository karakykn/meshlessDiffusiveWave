clc; clear all;
startX = datetime(2008,1,1,0,0,0);
Meshless_illinois_in = make_Meshless_df('../solution-ydk/segment2/run', 17, startX);
illinois_reg_in = regularize_to_hourly(Meshless_illinois_in);


% Define date range
startDate = datetime(2008,6,1,0,0,0);
endDate   = datetime(2008,9,1,0,0,0);
mask = illinois_reg_in.Time >= startDate & illinois_reg_in.Time <= endDate;
% Filter
writetable(illinois_reg_in(mask,:),'illinois_inflowbc.xlsx')

Meshless_miss_in = make_Meshless_df('../solution-ydk/segment1/run', 21, startX);
miss_reg_in = regularize_to_hourly(Meshless_miss_in);
% miss_reg_in = miss_reg_in(mask,:);
writetable(miss_reg_in(mask,:),'miss_inflowbc.xlsx')
%% Functions
function T = make_Meshless_df(run_dir, node_id, start_date)
    % Collect Q values
    [secs, dis] = collect_Q_values(run_dir, node_id);

    % Sort by seconds (safety, like Python code)
    [secs, order] = sort(secs);
    dis = dis(order);

    % Create table (similar to pandas DataFrame)
    T = table(secs, dis, ...
        'VariableNames', {'seconds', 'discharge_cms'});

    % Create Date column
    T.Date = start_date + seconds(T.seconds);
end

function [times, values] = collect_Q_values(base_dir, node, verbose)
%COLLECT_Q_VALUES Collect values from 'Q' files in numeric subdirectories
%
% Returns sorted arrays (times, values).

    if nargin < 1 || isempty(base_dir)
        base_dir = 'run';
    end
    if nargin < 2
        node = -1;
    end
    if nargin < 3
        verbose = false;
    end

    % MATLAB is 1-indexed; Python code assumes 0-indexed
    if node >= 0
        node = node + 1;
    end

    % Gather numeric subdirectories
    d = dir(base_dir);
    numeric_dirs = [];

    for i = 1:length(d)
        if ~d(i).isdir || ismember(d(i).name, {'.', '..'})
            continue
        end

        t = str2double(d(i).name);
        if ~isnan(t)
            numeric_dirs = [numeric_dirs; t, i]; %#ok<AGROW>
        else
            if verbose
                fprintf('Skipping non-numeric directory: %s\n', d(i).name);
            end
        end
    end

    % Sort by numeric time
    numeric_dirs = sortrows(numeric_dirs, 1);

    times   = [];
    values  = [];
    skipped = {};

    for k = 1:size(numeric_dirs, 1)
        t = numeric_dirs(k, 1);
        dir_idx = numeric_dirs(k, 2);
        name = d(dir_idx).name;

        qpath = fullfile(base_dir, name, 'Q');

        if ~isfile(qpath)
            skipped(end+1,:) = {name, 'no Q file'}; %#ok<AGROW>
            if verbose
                fprintf('Skipping %s: Q file not found\n', name);
            end
            continue
        end

        try
            data = load(qpath);
        catch ME
            skipped(end+1,:) = {name, ['load error: ' ME.message]}; %#ok<AGROW>
            if verbose
                fprintf('Skipping %s: error loading Q: %s\n', name, ME.message);
            end
            continue
        end

        if isempty(data)
            skipped(end+1,:) = {name, 'empty Q'}; %#ok<AGROW>
            if verbose
                fprintf('Skipping %s: Q is empty\n', name);
            end
            continue
        end

        try
            if isvector(data)
                val = data(node);
            else
                val = data(node, end);
            end
        catch ME
            skipped(end+1,:) = {name, ['indexing error: ' ME.message]}; %#ok<AGROW>
            if verbose
                fprintf('Skipping %s: indexing error: %s\n', name, ME.message);
            end
            continue
        end

        times(end+1,1)  = t;      %#ok<AGROW>
        values(end+1,1) = val;    %#ok<AGROW>
    end

    if verbose && ~isempty(skipped)
        fprintf('Skipped %d entries. Example skips:\n', size(skipped,1));
        disp(skipped(1:min(5,end), :))
    end
end

function T_out = regularize_to_hourly(T_in)

    T_in.Time = datetime(T_in{:,3}, 'InputFormat', 'dd-MMM-yyyy HH:mm:ss');
    TT = timetable(T_in.Time, T_in{:,2}, 'VariableNames', {'Flow'});
    
    newTime = (TT.Time(1):hours(1):TT.Time(end))';
    TT_hourly = retime(TT, newTime, 'linear');
    
    T_out = timetable2table(TT_hourly);

end