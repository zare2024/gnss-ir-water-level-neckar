%% GNSS-IR water level: Neckar river, Ober-Esslingen
% Reflector height (antenna height above the water surface) from L1 SNR
% observations of GPS, Galileo and GLONASS satellites.
%
% Each satellite arc over the river is detrended to isolate the multipath
% oscillation in SNR, its dominant frequency in sin(elevation) is found with
% a Lomb-Scargle periodogram, and the frequency is converted to a height.
% Arc estimates are combined with a median-based outlier window.
%
% Author: Shaghayegh Zare. Developed from Lab 3 of "Satellite Geodesy
% Observation Techniques" (University of Stuttgart, WS 2025/26, team work
% with Haibo and Shiwei); code rewritten and restructured for this repository.
%
% Requires MATLAB R2016b or newer with the Signal Processing Toolbox (plomb).
% Input: gnss_ir.mat next to this script (not included, see README).

clear; close all; clc;

%% Settings
cfg.dataFile      = 'gnss_ir.mat';
cfg.wavelength    = 0.19029367;   % GPS L1 carrier wavelength (m)
cfg.azimuthRange  = [270 330];    % deg, sector facing the river
cfg.heightRange   = [0.1 500];    % m, reflector heights searched
cfg.nFrequencies  = 2000;         % periodogram frequencies
cfg.polyOrder     = 2;            % trend model for the direct signal
cfg.minSamples    = 100;          % minimum samples per arc
cfg.outlierWindow = 1.5;          % m, accepted distance from the median
cfg.figureDir     = 'figures';

%% Estimate one reflector height per satellite arc
S = load(cfg.dataFile, 'gnss_ir');
tracks = estimateArcs(S.gnss_ir, cfg);
if height(tracks) == 0
    error('No satellite arc passed the quality checks.');
end

%% Combine arcs and report
result = combineArcs(tracks.Height_m, cfg.outlierWindow);
printReport(tracks, result, cfg);
plotArcs(tracks, result, cfg);


%% ------------------------------------------------------------------------
%  Local functions
%  ------------------------------------------------------------------------

function tracks = estimateArcs(data, cfg)
% One row per satellite arc that passes the quality checks.
    sats = fieldnames(data);
    satCol = {}; sysCol = {}; arcCol = {}; hCol = []; nCol = [];

    for k = 1:numel(sats)
        T = data.(sats{k});
        if ismember('time', T.Properties.VariableNames)
            t = T.time;
        else
            t = (1:height(T))';
        end

        inSector = T.azimut >= cfg.azimuthRange(1) & T.azimut <= cfg.azimuthRange(2);
        if nnz(inSector) < cfg.minSamples
            continue
        end

        [~, order] = sort(t(inSector));
        elev = T.elevation(inSector);  elev = elev(order);
        snr  = T.snr1(inSector);       snr  = snr(order);

        arcs = splitAtZenithPass(elev);
        for a = 1:size(arcs, 1)
            idx = arcs{a, 2};
            if numel(idx) < cfg.minSamples
                continue
            end
            [Hr, nUsed] = reflectorHeight(elev(idx), snr(idx), cfg);
            if isnan(Hr)
                continue
            end
            satCol{end+1, 1} = sats{k};              %#ok<AGROW>
            sysCol{end+1, 1} = systemName(sats{k}); %#ok<AGROW>
            arcCol{end+1, 1} = arcs{a, 1};           %#ok<AGROW>
            hCol(end+1, 1)   = Hr;                   %#ok<AGROW>
            nCol(end+1, 1)   = nUsed;                %#ok<AGROW>
        end
    end

    tracks = table(satCol, sysCol, arcCol, hCol, nCol, 'VariableNames', ...
        {'Satellite', 'System', 'Arc', 'Height_m', 'Samples'});
end


function arcs = splitAtZenithPass(elev)
% Rising and setting part of a pass, split at the highest elevation.
    [~, iTop] = max(elev);
    arcs = {'Ascending', 1:iTop; 'Descending', iTop:numel(elev)};
end


function [Hr, nUsed] = reflectorHeight(elevDeg, snr, cfg)
% Reflector height from the dominant SNR oscillation in sin(elevation).
% Returns NaN when too few distinct elevation samples remain.
    Hr = NaN;
    [x, order] = sort(sind(elevDeg));
    y = snr(order);
    [x, first] = unique(x);          % keep the first of repeated values
    y = y(first);
    nUsed = numel(x);
    if nUsed < cfg.minSamples
        return
    end

    trend = polyval(polyfit(x, y, cfg.polyOrder), x);   % direct signal
    fGrid = linspace(2 * cfg.heightRange(1) / cfg.wavelength, ...
                     2 * cfg.heightRange(2) / cfg.wavelength, cfg.nFrequencies);
    [power, f] = plomb(y - trend, x, fGrid);
    [~, iPeak] = max(power);
    Hr = f(iPeak) * cfg.wavelength / 2;                 % f = 2 Hr / lambda
end


function r = combineArcs(h, window)
% Median-centred window, then mean and standard deviation of the kept arcs.
    r.median  = median(h);
    r.kept    = abs(h - r.median) < window;
    r.mean    = mean(h(r.kept));
    r.std     = std(h(r.kept));
    r.window  = window;
end


function name = systemName(prn)
    switch upper(prn(1))
        case 'G', name = 'GPS';
        case 'E', name = 'Galileo';
        case 'R', name = 'GLONASS';
        case 'C', name = 'BeiDou';
        otherwise, name = 'other';
    end
end


function printReport(tracks, r, cfg)
    fprintf('Arcs estimated: %d (azimuth %d-%d deg)\n\n', height(tracks), ...
        cfg.azimuthRange(1), cfg.azimuthRange(2));
    fprintf('%-5s %-8s %-11s %10s %8s  %s\n', 'Sat', 'System', 'Arc', 'Hr (m)', 'Samples', 'Status');
    for i = 1:height(tracks)
        if r.kept(i), status = 'kept'; else, status = 'outlier'; end
        fprintf('%-5s %-8s %-11s %10.2f %8d  %s\n', tracks.Satellite{i}, tracks.System{i}, ...
            tracks.Arc{i}, tracks.Height_m(i), tracks.Samples(i), status);
    end

    fprintf('\nMedian of all arcs:   %.2f m\n', r.median);
    fprintf('Window:               +/- %.2f m\n', r.window);
    fprintf('Arcs kept / removed:  %d / %d\n', nnz(r.kept), nnz(~r.kept));

    fprintf('\nPer constellation (kept arcs):\n');
    systems = unique(tracks.System(r.kept));
    for s = 1:numel(systems)
        sel = r.kept & strcmp(tracks.System, systems{s});
        fprintf('  %-8s %d arcs, mean Hr = %.2f m\n', systems{s}, nnz(sel), mean(tracks.Height_m(sel)));
    end

    fprintf('\nREFLECTOR HEIGHT:     %.2f +/- %.2f m\n', r.mean, r.std);
end


function plotArcs(tracks, r, cfg)
% Kept arcs with the combined estimate and its standard deviation.
    kept = tracks(r.kept, :);
    x = 1:height(kept);
    labels = strcat(kept.Satellite, '-', kept.Arc);
    span = [0.5, numel(x) + 0.5];

    fig = figure('Color', 'w', 'Position', [100 100 760 380]);
    hold on
    fill([span, fliplr(span)], [r.mean - r.std, r.mean - r.std, r.mean + r.std, r.mean + r.std], ...
         [0.85 0.9 1], 'EdgeColor', 'none');
    plot(span, [r.mean r.mean], 'Color', [0.1 0.3 0.7], 'LineWidth', 1.5);
    plot(x, kept.Height_m, 'o', 'MarkerSize', 7, 'MarkerFaceColor', [0.95 0.5 0.1], ...
         'MarkerEdgeColor', 'k');
    hold off
    set(gca, 'XTick', x, 'XTickLabel', labels, 'XLim', span);
    xtickangle(45);
    ylabel('Reflector height H_r (m)');
    title(sprintf('Arcs kept: H_r = %.2f \\pm %.2f m (%d arcs)', r.mean, r.std, numel(x)));
    grid on; box on

    if ~exist(cfg.figureDir, 'dir')
        mkdir(cfg.figureDir);
    end
    saveas(fig, fullfile(cfg.figureDir, 'arc_heights.png'));
end
