clear; clc; close all;

fs = 44100;
N = 1000;
sounds = {'e', 'a', 's', 't', 'x'};
fileNames = strcat(sounds, '.wav');

rootDir = fileparts(mfilename('fullpath'));
audioDir = fullfile(rootDir, 'audio');
if ~exist(audioDir, 'dir')
    audioDir = rootDir;
end

% Predictor orders and adaptation gains to test
pVals = [2, 4, 8, 16, 32];
muVals = [0.01, 0.05, 0.1, 0.3];

ssWin = 200;  % steady-state window

fprintf('Orders p: '); fprintf('%d ', pVals); fprintf('\n');
fprintf('Mu values: '); fprintf('%.4f ', muVals); fprintf('\n\n');

results = struct();

for si = 1:numel(fileNames)
    filePath = fullfile(audioDir, fileNames{si});
    [x, ok] = load_audio_segment(filePath, fs, N);
    if ~ok
        fprintf('Missing or unreadable file: %s (skipping)\n', filePath);
        continue;
    end

    fprintf('\nSound "%s":\n', sounds{si});

    best = struct('p', NaN, 'mu', NaN, 'mse', Inf, 'predGain', -Inf);
    table = [];

    for p = pVals
        for mu = muVals
            [~, e] = lms_predictor_fixed(x, p, mu);
            mse = mean(e(end-ssWin+1:end) .^ 2);
            predGain = 10 * log10(var(x(p+1:end), 1) / max(mse, eps));
            table = [table; p, mu, mse, predGain]; %#ok<AGROW>

            if mse < best.mse
                best.p = p;
                best.mu = mu;
                best.mse = mse;
                best.predGain = predGain;
            end
        end
    end

    % Gear shifting test (use best p from fixed-mu sweep)
    Px = mean(x .^ 2);
    cfg.muMax = 0.2;
    cfg.muMin = 0.005;
    cfg.alphaErr = 0.98;
    cfg.betaMu = 0.90;
    cfg.errLow = 0.02 * var(x, 1);
    cfg.errHigh = 0.20 * var(x, 1);

    [~, eGear, ~, muHist] = lms_predictor_gear(x, best.p, cfg);
    mseGear = mean(eGear(end-ssWin+1:end) .^ 2);
    predGainGear = 10 * log10(var(x(best.p+1:end), 1) / max(mseGear, eps));

    fprintf('Best fixed-mu: p=%d, mu=%.4f, SS-MSE=%.6f, predGain=%.2f dB\n', ...
        best.p, best.mu, best.mse, best.predGain);
    fprintf('Gear shift:    p=%d, mu in [%.4f, %.4f], SS-MSE=%.6f, predGain=%.2f dB\n', ...
        best.p, cfg.muMin, cfg.muMax, mseGear, predGainGear);

    results.(sounds{si}).table = table;
    results.(sounds{si}).best = best;
    results.(sounds{si}).gear = struct('mse', mseGear, 'predGain', predGainGear, ...
        'muMin', cfg.muMin, 'muMax', cfg.muMax, 'muHist', muHist);
end

% Error-curve plots for all sounds
figure('Name', 'Error curves for different \mu (all sounds)');
tl = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
pPlot = pVals(1);
firstLegend = true;

for si = 1:numel(sounds)
    if ~isfield(results, sounds{si})
        continue;
    end
    filePath = fullfile(audioDir, [sounds{si} '.wav']);
    [x, ok] = load_audio_segment(filePath, fs, N);
    if ~ok
        continue;
    end

    nexttile;
    hold on;
    for mu = muVals
        [~, e] = lms_predictor_fixed(x, pPlot, mu);
        plot(movmean(e .^ 2, 50), 'LineWidth', 1.1, ...
            'DisplayName', sprintf('\\mu=%.3f', mu));
    end
    hold off;
    grid on;
    xlabel('n', 'FontSize', 12);
    ylabel('Smoothed e^2[n]', 'FontSize', 12);
    title(sprintf('Sound "%s" (p=%d)', sounds{si}, pPlot), 'FontSize', 12);
    set(gca, 'FontSize', 11);
    if firstLegend
        legend('Location', 'best');
        firstLegend = false;
    end
end
sgtitle(tl, 'Prediction error curves for all sounds', 'FontSize', 14);

% Time-series plots for all sounds
figure('Name', 'Time Series (all sounds)');
tl2 = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for si = 1:numel(sounds)
    filePath = fullfile(audioDir, [sounds{si} '.wav']);
    [x, ok] = load_audio_segment(filePath, fs, N);
    if ~ok
        continue;
    end
    nexttile;
    plot((0:N-1) / fs, x, 'LineWidth', 1.0);
    grid on;
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('Amplitude', 'FontSize', 12);
    title(sprintf('Sound "%s"', sounds{si}), 'FontSize', 12);
    set(gca, 'FontSize', 11);
end
sgtitle(tl2, 'Time-series segments used for prediction (N=1000)', 'FontSize', 14);

% Gear-shift mu trajectories for all sounds, with reference mu values
muTypical = median(muVals);
figure('Name', 'Gear-Shift \mu Trajectories (all sounds)');
tlMu = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for si = 1:numel(sounds)
    if ~isfield(results, sounds{si})
        continue;
    end
    muHist = results.(sounds{si}).gear.muHist;
    muMin = results.(sounds{si}).gear.muMin;
    muMax = results.(sounds{si}).gear.muMax;

    nexttile;
    plot(muHist, 'LineWidth', 1.1);
    hold on;
    for mv = muVals
        yline(mv, ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    end
    yline(muTypical, '--', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.0, ...
        'DisplayName', sprintf('median(\\mu)=%.3f', muTypical));
    yline(muMin, '--', 'Color', [0.2 0.6 0.2], 'HandleVisibility', 'off');
    yline(muMax, '--', 'Color', [0.8 0.2 0.2], 'HandleVisibility', 'off');
    hold off;
    grid on;
    xlabel('n', 'FontSize', 12);
    ylabel('\\mu[n]', 'FontSize', 12);
    title(sprintf('Sound "%s"', sounds{si}), 'FontSize', 12);
    set(gca, 'FontSize', 11);
end
sgtitle(tlMu, 'Gear-shifting adaptation gain \\mu[n] (with muVals references)', 'FontSize', 14);

maxLag = 25;
for si = 1:numel(sounds)
    filePath = fullfile(audioDir, [sounds{si} '.wav']);
    [x, ok] = load_audio_segment(filePath, fs, N);
    if ~ok
        continue;
    end

    [pacfVals, lags] = pacf_ols(x, maxLag);

    figure('Name', sprintf('PACF: %s', sounds{si}));
    stem(lags, pacfVals, 'filled');
    grid on;
    xlabel('Lag', 'FontSize', 12);
    ylabel('PACF', 'FontSize', 12);
    title(sprintf('Partial Autocorrelation: "%s"', sounds{si}), 'FontSize', 12);
    set(gca, 'FontSize', 11);
end

pFixed = pVals(1);
muFixed = muVals(1);

figMu = figure('Name', 'Coefficient Convergence (varying \mu)');
tlMu = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
colorsMu = lines(numel(muVals));
firstLegendMu = true;

for si = 1:numel(sounds)
    filePath = fullfile(audioDir, [sounds{si} '.wav']);
    [x, ok] = load_audio_segment(filePath, fs, N);
    if ~ok
        continue;
    end

    ax = nexttile(tlMu);
    colororder(ax, colorsMu);
    hold(ax, 'on');
    for mi = 1:numel(muVals)
        mu = muVals(mi);
        [~, ~, W] = lms_predictor_fixed(x, pFixed, mu);
        c = colorsMu(mi, :);
        plot(ax, W(1, :), '-', 'LineWidth', 1.0, 'Color', c, ...
            'DisplayName', sprintf('\\mu=%.3f', mu));
        plot(ax, W(2, :), '--', 'LineWidth', 1.0, 'Color', c, ...
            'HandleVisibility', 'off');
    end
    hold(ax, 'off');
    grid(ax, 'on');
    xlabel(ax, 'n', 'FontSize', 14);
    ylabel(ax, 'Coefficient value', 'FontSize', 14);
    title(ax, sprintf('Sound "%s" (p=%d)', sounds{si}, pFixed), 'FontSize', 14);
    set(ax, 'FontSize', 12);
    if firstLegendMu
        legend(ax, 'Location', 'best');
        firstLegendMu = false;
    end
end
sgtitle(tlMu, 'Coefficient convergence (varying \mu) — solid a_1, dashed a_2', 'FontSize', 14);

figP = figure('Name', 'Coefficient Convergence (varying p)');
tlP = tiledlayout(3, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
colorsP = lines(numel(pVals));
firstLegendP = true;

for si = 1:numel(sounds)
    filePath = fullfile(audioDir, [sounds{si} '.wav']);
    [x, ok] = load_audio_segment(filePath, fs, N);
    if ~ok
        continue;
    end

    ax = nexttile(tlP);
    colororder(ax, colorsP);
    hold(ax, 'on');
    for pi = 1:numel(pVals)
        p = pVals(pi);
        [~, ~, W] = lms_predictor_fixed(x, p, muFixed);
        c = colorsP(pi, :);
        plot(ax, W(1, :), '-', 'LineWidth', 1.0, 'Color', c, ...
            'DisplayName', sprintf('p=%d', p));
        plot(ax, W(2, :), '--', 'LineWidth', 1.0, 'Color', c, ...
            'HandleVisibility', 'off');
    end
    hold(ax, 'off');
    grid(ax, 'on');
    xlabel(ax, 'n', 'FontSize', 14);
    ylabel(ax, 'Coefficient value', 'FontSize', 14);
    title(ax, sprintf('Sound "%s" (\\mu=%.3f)', sounds{si}, muFixed), 'FontSize', 14);
    set(ax, 'FontSize', 12);
    if firstLegendP
        legend(ax, 'Location', 'best');
        firstLegendP = false;
    end
end
sgtitle(tlP, 'Coefficient convergence (varying p) — solid a_1, dashed a_2', 'FontSize', 14);

function [yhat, e, W] = lms_predictor_fixed(x, p, mu)
    N = numel(x);
    yhat = zeros(N, 1);
    e = zeros(N, 1);
    W = zeros(p, N);
    w = zeros(p, 1);

    for n = p+1:N
        xvec = x(n-1:-1:n-p);
        yhat(n) = w.' * xvec;
        e(n) = x(n) - yhat(n);
        w = w + mu * e(n) * xvec;
        W(:, n) = w;
    end
end

function [yhat, e, W, muHist] = lms_predictor_gear(x, p, cfg)
    N = numel(x);
    yhat = zeros(N, 1);
    e = zeros(N, 1);
    W = zeros(p, N);
    muHist = zeros(N, 1);
    w = zeros(p, 1);
    muNow = cfg.muMax;
    errPow = cfg.errHigh;

    for n = p+1:N
        xvec = x(n-1:-1:n-p);
        yhat(n) = w.' * xvec;
        e(n) = x(n) - yhat(n);

        errPow = cfg.alphaErr * errPow + (1 - cfg.alphaErr) * (e(n) ^ 2);
        if errPow <= cfg.errLow
            muTarget = cfg.muMin;
        elseif errPow >= cfg.errHigh
            muTarget = cfg.muMax;
        else
            ratio = (errPow - cfg.errLow) / (cfg.errHigh - cfg.errLow);
            muTarget = cfg.muMin + ratio * (cfg.muMax - cfg.muMin);
        end
        muNow = cfg.betaMu * muNow + (1 - cfg.betaMu) * muTarget;
        muNow = min(max(muNow, cfg.muMin), cfg.muMax);

        w = w + muNow * e(n) * xvec;
        W(:, n) = w;
        muHist(n) = muNow;
    end
end

function [x, ok] = load_audio_segment(filePath, fs, N)
    ok = false;
    if ~exist(filePath, 'file')
        x = [];
        return;
    end

    [x, fsIn] = audioread(filePath);
    if size(x, 2) > 1
        x = mean(x, 2);
    end
    if fsIn ~= fs
        if exist('resample', 'file') == 2
            x = resample(x, fs, fsIn);
        else
            tOld = (0:numel(x)-1) / fsIn;
            tNew = (0:round(numel(x) * fs / fsIn) - 1) / fs;
            x = interp1(tOld, x, tNew, 'linear', 'extrap');
        end
    end

    if numel(x) < N
        x = [x; zeros(N - numel(x), 1)];
    elseif numel(x) > N
        win = N;
        energy = conv(x.^2, ones(win, 1), 'valid');
        [~, idx] = max(energy);
        x = x(idx:idx+N-1);
    end
    x = x - mean(x);
    ok = true;
end

function [pacfVals, lags] = pacf_ols(x, maxLag)
    % Estimate PACF via Yule-Walker (OLS on autocorrelation)
    x = x(:);
    N = numel(x);
    r = xcorr(x, maxLag, 'biased');
    r = r(maxLag+1:end); % lags 0..maxLag

    pacfVals = zeros(maxLag, 1);
    for k = 1:maxLag
        R = toeplitz(r(1:k));
        rho = r(2:k+1);
        a = R \ rho;
        pacfVals(k) = a(end);
    end
    lags = 1:maxLag;
end
