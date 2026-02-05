% Xue R. et al., "Force loading on molecular clutches governs the stability of cell lamellipodia".
%
% Description:
%   This script performs a comprehensive analysis of oscillation data stored
%   in an Excel file. It is divided into two main parts:
%
%   Part 1: Statistical Analysis
%       - Extracts oscillation amplitude and period.
%       - Generates bar plots with individual data points (scatter).
%       - Performs statistical significance testing (Mann-Whitney U test).
%
%   Part 2: Spectral and Temporal Analysis
%       - Computes Power Spectral Density (PSD).
%       - Generates aggregate PSD plots with standard error shading.
%
% Input:
%   - 'All_Data_Combined.xlsx': Excel file containing time-series data.
%     Each sheet represents a sample. Sheet names must contain specific keywords.
%
% Output:
%   - Figures for statistical comparisons, representative traces, and PSDs.
%
% Date: Feb 2026

    clc; clear; close all;

    %% 1. Initialization and Global Settings
    filename = 'All_Data_Combined.xlsx'; 
    
    % Check for data file existence
    if ~isfile(filename)
        error('Error: Data file "%s" not found in the current directory.', filename); 
    end
    
    % Define Group Keywords (Used to parse Excel sheet names)
    groupKeywords = {'MIIB-KO_MIIA-KD', 'Blebb', 'MIIB-KO', 'Control', 'Calyculin A'};
    
    % Define Group Labels (Used for Legends/Figures)
    groupLabels = {'MIIB-KO/MIIA-KD', 'Blebbistatin', 'MIIB-KO', 'Control', 'Calyculin A'};

    fprintf('=== Starting Comprehensive Analysis ===\n');
    fprintf('Data Source: %s\n\n', filename);

    %% ====================================================================
    %% PART 1: Statistical Bar Plots (Amplitude & Period)
    %% ====================================================================
    fprintf('--- Part 1: Generating Statistical Bar Plots ---\n');
    
    nGroups = length(groupKeywords);
    
    % 1.1 Color Configuration (Sequential Blues for Bar Plots)
    % Attempt to use 'othercolor' if installed, otherwise use a simulation function.
    try
        barColors = othercolor('Blues7', 10 * 5);
    catch
        barColors = generate_blues7_sim(10 * nGroups);
    end
    
    % 1.2 Data Extraction
    stats_Amp = cell(1, nGroups); 
    stats_Period = cell(1, nGroups); 
    
    [~, sheetList] = xlsfinfo(filename);
    
    for g = 1:nGroups
        key = groupKeywords{g};
        
        % Filter sheets based on keywords
        idx = contains(sheetList, key, 'IgnoreCase', true);
        
        % Exclusion logic for overlapping naming conventions
        if strcmp(key, 'MIIB-KO'), idx = idx & ~contains(sheetList, 'KD', 'IgnoreCase', true); end
        if strcmp(key, 'Blebb'),   idx = idx & ~contains(sheetList, 'BBI', 'IgnoreCase', true); end
        targetSheets = sheetList(idx);
        
        raw_Amps = []; 
        raw_Periods = [];
        
        for k = 1:length(targetSheets)
            try
                % Read and preprocess data
                raw = readmatrix(filename, 'Sheet', targetSheets{k});
                raw = rmmissing(raw);
                if size(raw, 2) < 2, continue; end
                
                time_vec = raw(:, 1); 
                y = raw(:, 2);
                
                % Sort and remove duplicates
                [time_vec, I] = sort(time_vec); y = y(I); 
                [time_vec, I] = unique(time_vec); y = y(I);
                if length(time_vec) < 10, continue; end
                
                % Detrending (Polynomial fit order 3)
                p = polyfit(time_vec, y, 3);
                y_detrend = y - polyval(p, time_vec);
                
                % Dynamic Thresholding Strategy
                % If signal noise/fluctuation is high, use higher threshold (0.1).
                % Otherwise, detect smaller peaks (0.0).
                sig_strength = std(y_detrend); 
                if sig_strength > 0.1 
                    current_thresh = 0.1;
                else
                    current_thresh = 0.0;
                end
                
                % Extract features
                [amps, pers, ~, ~, ~, ~] = extract_features(time_vec, y_detrend, current_thresh);
                
                raw_Amps = [raw_Amps, amps];
                raw_Periods = [raw_Periods, pers(:)'];
            catch
                % Skip sheet on read error
                continue; 
            end
        end
        
        % Remove outliers (using Quartile method) and store
        if ~isempty(raw_Amps)
            stats_Amp{g} = raw_Amps(~isoutlier(raw_Amps, 'quartiles')); 
        end
        if ~isempty(raw_Periods)
            stats_Period{g} = raw_Periods(~isoutlier(raw_Periods, 'quartiles')); 
        end
    end
    
    % 1.3 Visualization (Bar Plots)
    % Amplitude: Y-limit [0, 4]
    plot_custom_bar('Fig1_Amp_Stats', stats_Amp, barColors, [0, 4]); 
    % Period: Y-limit [0, 80]
    plot_custom_bar('Fig2_Period_Stats', stats_Period, barColors, [0, 80]);
    
    fprintf('Part 1 Complete. Statistical figures generated.\n\n');


    %% ====================================================================
    %% PART 2: PSD Analysis & Representative Traces
    %% ====================================================================
    fprintf('--- Part 2: Generating PSD and Representative Traces ---\n');

    % 2.1 Parameters
    poly_order = 3;             % Order of polynomial for detrending
    target_Fs = 10;             % Resampling frequency (Hz)
    freq_axis = 0:0.002:0.1;    % Frequency axis for PSD
    
    % 2.2 Color Configuration (Categorical Colors for PSD)
    % Distinct colors for clear legend differentiation
    psdColors = [
        0.85, 0.65, 0.10;   % Mustard Gold (MIIB-KO/MIIA-KD)
        0.05, 0.35, 0.75;   % Deep Blue (Blebbistatin)
        0.85, 0.25, 0.10;   % Red-Orange (MIIB-KO)
        0.40, 0.70, 0.20;   % Green (Control)
        0.55, 0.20, 0.65    % Purple (Calyculin A)
    ];

    % 2.3 Container Initialization
    stats_PSD_All = cell(1, nGroups);
    
    % 2.4 Main Loop (Trace Extraction & PSD Calculation)
    for g = 1:nGroups
        key = groupKeywords{g};
        thisColor = psdColors(g, :);
        
        % Filter Sheets (Same logic as Part 1)
        idx = contains(sheetList, key, 'IgnoreCase', true);
        if strcmp(key, 'MIIB-KO'), idx = idx & ~contains(sheetList, 'KD', 'IgnoreCase', true); end
        if strcmp(key, 'Blebb'),   idx = idx & ~contains(sheetList, 'BBI', 'IgnoreCase', true); end
        targetSheets = sheetList(idx);
        
        if isempty(targetSheets), continue; end
        
        group_PSDs = [];
        plotBuffer = struct([]); % Buffer to store first 3 samples for plotting
        count = 0; 
        
        for k = 1:length(targetSheets)
            try
                % Read and Preprocess
                raw = readmatrix(filename, 'Sheet', targetSheets{k});
                raw = rmmissing(raw);
                if size(raw, 2) < 2, continue; end
                
                time_vec = raw(:, 1); y = raw(:, 2);
                [time_vec, I] = sort(time_vec); y = y(I); 
                [time_vec, I] = unique(time_vec); y = y(I);
                if length(time_vec) < 10, continue; end
                
                % Normalize time and amplitude start points
                time_vec = time_vec - time_vec(1); 
                y = y - y(1); 
                
                % Detrend
                p = polyfit(time_vec, y, poly_order);
                trend = polyval(p, time_vec);
                y_detrend = y - trend;
                
                % Calculate PSD (with resampling)
                t_u = time_vec(1):(1/target_Fs):time_vec(end);
                y_u = interp1(time_vec, y, t_u, 'linear');
                p_u = polyfit(t_u, y_u, poly_order);
                y_u_dt = y_u - polyval(p_u, t_u);
                
                [pxx, ~] = periodogram(y_u_dt, rectwin(length(y_u_dt)), freq_axis, target_Fs);
                group_PSDs = [group_PSDs, pxx(:)];
                
                % Cache first 3 valid samples for Trace Plotting
                if count < 3
                    run_thresh = 0.1;
                    [~, ~, pks, locs, vlys, vlylocs] = extract_features(time_vec, y_detrend, run_thresh);
                    % Fallback: if no peaks found, lower threshold for visualization
                    if length(pks) < 3
                         [~, ~, pks, locs, vlys, vlylocs] = extract_features(time_vec, y_detrend, 0.0);
                    end
                    
                    count = count + 1;
                    plotBuffer(count).t = time_vec;
                    plotBuffer(count).y = y;
                    plotBuffer(count).trend = trend;
                    plotBuffer(count).y_detrend = y_detrend;
                    plotBuffer(count).pks = pks;
                    plotBuffer(count).locs = locs;
                    plotBuffer(count).vlys = vlys;
                    plotBuffer(count).vlylocs = vlylocs;
                end
            catch
                continue;
            end
        end
        stats_PSD_All{g} = group_PSDs;
        
        % ---------------------------------------------------------
        % Generate Individual Trace Plots (Raw + Detrended)
        % ---------------------------------------------------------
        if isempty(plotBuffer), continue; end
        
        % Calculate unified Y-limits for this group for consistency
        min_raw_val = inf; max_raw_val = -inf;
        max_detrend_abs = 0;
        
        for i = 1:length(plotBuffer)
            min_raw_val = min(min_raw_val, min(plotBuffer(i).y));
            max_raw_val = max(max_raw_val, max(plotBuffer(i).y));
            max_detrend_abs = max(max_detrend_abs, max(abs(plotBuffer(i).y_detrend)));
        end
        
        ylim_raw = [floor(min_raw_val), ceil(max_raw_val)]; 
        ceil_dt = ceil(max_detrend_abs);
        if ceil_dt == 0, ceil_dt = 1; end
        ylim_detrend = [-ceil_dt, ceil_dt];
        
        for i = 1:length(plotBuffer)
            pb = plotBuffer(i);
            
            % Round X-axis limit to nearest 50
            xlim_sample = ceil(max(pb.t) / 50) * 50;
            if xlim_sample == 0, xlim_sample = 50; end 
            
            f_name = sprintf('Trace_Group%d_%s_Sample%d', g, key, i);
            figure('Name', f_name, 'Color', 'w', 'Position', [500, 500, 190, 200]);
            t = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');
            
            % Top Plot: Raw Data + Trend
            ax1 = nexttile;
            plot(pb.t, pb.y, 'Color', [0.1 0.1 0.1], 'LineWidth', 1); hold on;
            plot(pb.t, pb.trend, 'r--', 'LineWidth', 1.5);
            xlim([0, xlim_sample]); ylim(ylim_raw);
            box on; set(ax1, 'XTickLabel', {}, 'FontName', 'Arial', 'FontSize', 10);
            
            % Bottom Plot: Detrended Data + Detected Peaks/Valleys
            ax2 = nexttile;
            plot(pb.t, pb.y_detrend, 'Color', thisColor, 'LineWidth', 1.2); hold on;
            if ~isempty(pb.locs), scatter(pb.locs, pb.pks, 15, 'b', 'filled'); end
            if ~isempty(pb.vlylocs), scatter(pb.vlylocs, -pb.vlys, 15, 'g', 'filled'); end
            xlim([0, xlim_sample]); ylim(ylim_detrend);
            box on; set(ax2, 'FontName', 'Arial', 'FontSize', 10);
        end
    end
    
    % 2.5 Generate Aggregate PSD Plots
    % View 1: Zoomed out
    plot_clean_psd(freq_axis, stats_PSD_All, psdColors, [0, 0.06], [0, 14], 'Fig3_PSD_ZoomOut', groupLabels);
    % View 2: Zoomed in (Low Frequency)
    plot_clean_psd(freq_axis, stats_PSD_All, psdColors, [0.02, 0.06], [0, 3], 'Fig4_PSD_ZoomIn', groupLabels);
    
    fprintf('Part 2 Complete. Traces and PSD figures generated.\n');
    fprintf('=== Analysis Finished Successfully ===\n');



%% 1. Plotting Function: Custom Bar Plot (Part 1)
function plot_custom_bar(figName, dataCell, colors, y_limits)
    f = figure('Name', figName, 'Color', 'w', 'Position', [300, 300, 350, 320]);
    hold on;
    nGroups = length(dataCell);
    
    % Draw Bars and Scatter Points
    for i = 1:nGroups
        vals = dataCell{i};
        if isempty(vals), vals = 0; end
        m = mean(vals); 
        sd = std(vals);
        
        % Get color based on group index (sequential blue)
        curr_color = colors(10*i-4, :); 

        % Bar settings
        b = bar(i, m, 0.5);
        b.FaceColor = 'flat';       
        b.CData = curr_color;       
        b.EdgeColor = 'none';       
        
        % Scatter settings (Individual data points)
        jitter = (rand(size(vals))-0.5)*0.25;
        scatter(i+jitter, vals, "filled", "MarkerFaceColor", curr_color, 'SizeData', 8);

        % Error bar settings (Standard Deviation)
        errorbar(i, m, sd, "LineStyle", "none", 'Color', [0.15 0.15 0.15], 'LineWidth', 0.5);
    end
    
    % Axis adjustments
    xlim([0.5 5.5]);
    ylim(y_limits);
    set(gca, 'FontSize', 8, 'Color', 'none', 'TickLength', [0.02 0.02], 'Layer', 'top', 'TickDir', 'out');
    set(gcf, "Units", "centimeters", "Position", [30, 10, 5.4, 3.5]);
    xlabel(''); ylabel(''); title('');

    % Add Statistical Significance Lines
    y_max = y_limits(2);
    h_Top = y_max * 0.93;     
    h_Mid = y_max * 0.82;     
    h_Bot = y_max * 0.71;     

    % Level 1: Top comparisons
    add_sig_line(dataCell, 1, 3, h_Top);
    add_sig_line(dataCell, 3, 5, h_Top);

    % Level 2: Middle comparisons
    add_sig_line(dataCell, 2, 3, h_Mid);
    add_sig_line(dataCell, 3, 4, h_Mid);

    % Level 3: Bottom comparisons
    add_sig_line(dataCell, 1, 2, h_Bot); 
    add_sig_line(dataCell, 4, 5, h_Bot); 
end

%% 2. Statistical Testing: Mann-Whitney U & Significance Lines
function add_sig_line(dataCell, idx1, idx2, y_h)
    d1 = dataCell{idx1}; d2 = dataCell{idx2};
    if isempty(d1) || isempty(d2), return; end
    
    % Perform Mann-Whitney U test (Wilcoxon rank sum)
    p = ranksum(d1, d2);
    
    % Determine significance label
    if p > 0.05
        txt = 'NS'; fontSize = 8;
    else
        if p < 0.0001
            txt = '****';
        elseif p < 0.001
            txt = '***';
        elseif p < 0.01
            txt = '**';
        else
            txt = '*';
        end
        fontSize = 12;
    end
    
    % Draw bracket
    tick_len = y_h * 0.03; 
    x_vec = [idx1, idx1, idx2, idx2];
    y_vec = [y_h - tick_len, y_h, y_h, y_h - tick_len];
    
    plot(x_vec, y_vec, 'k-', 'LineWidth', 1);
    
    % Add text
    text((idx1 + idx2)/2, y_h + y_h*0.01, txt, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', fontSize, 'FontName', 'Arial', ...
        'FontWeight', 'bold', 'Color', 'k');
end

%% 3. Plotting Function: PSD (Part 2)
function plot_clean_psd(freq_axis, psd_data_cell, colors, x_lims, y_lims, figName, labels)
    figure('Name', figName, 'Color', 'w', 'Position', [300, 300, 350, 250]);
    hold on;
    
    plotted_handles = []; 
    plotted_labels = {};  
    
    for i = 1:length(psd_data_cell)
        psd_mat = psd_data_cell{i};
        if isempty(psd_mat), continue; end
        
        % Calculate Statistics
        m = mean(psd_mat, 2)';
        s = std(psd_mat, 0, 2)';
        sem = s/sqrt(size(psd_mat, 2));
        color = colors(i, :);
        
        % Plot Standard Error Shading (not in legend)
        fill([freq_axis fliplr(freq_axis)], [m+sem fliplr(m-sem)], ...
            color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
        % Plot Mean Line
        h = plot(freq_axis, m, 'Color', color, 'LineWidth', 2);
        
        plotted_handles = [plotted_handles, h];
        plotted_labels{end+1} = labels{i};
    end
    
    xlim(x_lims); ylim(y_lims);
    
    % Add Legend
    if ~isempty(plotted_handles)
        legend(plotted_handles, plotted_labels, ...
            'Box', 'off', 'Location', 'northeast', ...
            'FontName', 'Arial', 'FontSize', 10);
    end
    
    % Styling
    set(gca, 'Box', 'off', 'LineWidth', 1.2, 'FontSize', 10, 'FontName', 'Arial', 'TickDir', 'out');
    xlabel(''); ylabel(''); title('');
end

%% 4. Feature Extraction: Peak & Valley Detection
function [amps, pers, pks, locs, vlys, vlylocs] = extract_features(t, y, thresh)
    % Find Peaks (Local Maxima)
    [pks, locs] = findpeaks(y, t, 'MinPeakProminence', thresh);
    
    % Find Valleys (Local Minima via inversion)
    [vlys, vlylocs] = findpeaks(-y, t, 'MinPeakProminence', thresh);
    
    amps = []; pers = [];
    
    % Pair Peaks with subsequent Valleys to calculate amplitude
    for j = 1:length(locs)
        idx_v = find(vlylocs > locs(j), 1);
        if ~isempty(idx_v)
            real_amp = pks(j) + vlys(idx_v); % Total Amplitude (Peak + Valley depth)
            if real_amp > thresh, amps(end+1) = real_amp; end
        end
    end
    
    % Calculate Period (Time difference between consecutive peaks/valleys)
    if ~isempty(amps)
        pers = [diff(locs); diff(vlylocs)]'; 
    end
end

%% 5. Utility: Color Map Simulation
function cmap = generate_blues7_sim(n)
    % Fallback function if 'othercolor' library is not available.
    % Simulates a sequential Blue gradient.
    c_start = [0.8, 0.9, 1.0]; % Light Blue
    c_end   = [0.0, 0.2, 0.5]; % Dark Blue
    
    r = linspace(c_start(1), c_end(1), n)';
    g = linspace(c_start(2), c_end(2), n)';
    b = linspace(c_start(3), c_end(3), n)';
    cmap = [r, g, b];
end