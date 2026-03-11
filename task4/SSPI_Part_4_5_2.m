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

pVals = [2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32];
maxLag = 25;

figure('Name', 'PACF (all sounds)');
hold on;
colors = lines(numel(sounds));
for si = 1:numel(sounds)
    filePath = fullfile(audioDir, fileNames{si});
    [x, ok] = load_audio_segment(filePath, fs, N);
    if ~ok
        continue;
    end
    [pacfVals, lags] = pacf_ols(x, maxLag);
    h = stem(lags, pacfVals, 'filled', 'DisplayName', sprintf('"%s"', sounds{si}));
    set(h, 'Color', colors(si, :), 'MarkerFaceColor', colors(si, :));
end
hold off;
grid on;
xlabel('Lag', 'FontSize', 12);
ylabel('PACF', 'FontSize', 12);
title('Partial Autocorrelation (all sounds)', 'FontSize', 12);
legend('Location', 'best');

figure('Name', 'AIC and MDL vs Model Order (all sounds)');
hold on;
colors = lines(numel(sounds));
hList = gobjects(0);
lblList = {};
for si = 1:numel(sounds)
    filePath = fullfile(audioDir, fileNames{si});
    [x, ok] = load_audio_segment(filePath, fs, N);
    if ~ok
        continue;
    end
    [aicVals, mdlVals] = aic_mdl_ar(x, pVals);
    [~, idxA] = min(aicVals);
    [~, idxM] = min(mdlVals);
    fprintf('Sound "%s": AIC optimal p = %d, MDL optimal p = %d\n', ...
        sounds{si}, pVals(idxA), pVals(idxM));
    c = colors(si, :);
    h1 = plot(pVals, aicVals, 'o-', 'LineWidth', 1.2, 'Color', c);
    h2 = plot(pVals, mdlVals, 's--', 'LineWidth', 1.2, 'Color', c);
    hList(end+1) = h1; %#ok<AGROW>
    lblList{end+1} = sprintf('"%s" AIC', sounds{si}); %#ok<AGROW>
    hList(end+1) = h2; %#ok<AGROW>
    lblList{end+1} = sprintf('"%s" MDL', sounds{si}); %#ok<AGROW>
end
hold off;
grid on;
xlabel('Model order p', 'FontSize', 12);
ylabel('Criterion value', 'FontSize', 12);
title('AIC and MDL criteria vs model order (all sounds)', 'FontSize', 12);
legend(hList, lblList, 'Location', 'best');

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
    x = x(:);
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

function [aicVals, mdlVals] = aic_mdl_ar(x, pList)
    x = x(:);
    N = numel(x);
    pList = pList(:).';
    aicVals = zeros(size(pList));
    mdlVals = zeros(size(pList));
    for i = 1:numel(pList)
        p = pList(i);
        r = xcorr(x, p, 'biased');
        r = r(p+1:end); % lags 0..p
        R = toeplitz(r(1:p));
        rho = r(2:p+1);
        a = R \ rho;
        sigma2 = max(r(1) - rho.' * a, eps);
        aicVals(i) = log(sigma2) + (2 * p) / N;
        mdlVals(i) = log(sigma2) + (p * log(N)) / N;
    end
end
