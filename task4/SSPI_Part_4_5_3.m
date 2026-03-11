clear; clc; close all;

fs = 16000;
N = 1000;
sounds = {'e', 'a', 's', 't', 'x'};
fileNames = strcat(sounds, '.wav');

rootDir = fileparts(mfilename('fullpath'));
audioDir = fullfile(rootDir, 'audio');
if ~exist(audioDir, 'dir')
    audioDir = rootDir;
end

pVals = [2, 4, 8, 16];
muVals = [0.01, 0.05, 0.1, 0.3];
ssWin = 200;

fprintf('Prediction gain at fs=%d Hz, N=%d\n', fs, N);
fprintf('Orders p: '); fprintf('%d ', pVals); fprintf('\n');
fprintf('Mu values: '); fprintf('%.4f ', muVals); fprintf('\n\n');

for si = 1:numel(fileNames)
    filePath = fullfile(audioDir, fileNames{si});
    [x, ok] = load_audio_segment(filePath, fs, N);
    if ~ok
        fprintf('Missing or unreadable file: %s (skipping)\n', filePath);
        continue;
    end

    best = struct('p', NaN, 'mu', NaN, 'mse', Inf, 'predGain', -Inf);
    for p = pVals
        for mu = muVals
            [~, e] = lms_predictor_fixed(x, p, mu);
            mse = mean(e(end-ssWin+1:end) .^ 2);
            predGain = 10 * log10(var(x(p+1:end), 1) / max(mse, eps));
            if mse < best.mse
                best.p = p;
                best.mu = mu;
                best.mse = mse;
                best.predGain = predGain;
            end
        end
    end

    fprintf('Sound "%s": best p=%d, mu=%.3f, SS-MSE=%.6f, Rp=%.2f dB\n', ...
        sounds{si}, best.p, best.mu, best.mse, best.predGain);
end

function [yhat, e] = lms_predictor_fixed(x, p, mu)
    N = numel(x);
    yhat = zeros(N, 1);
    e = zeros(N, 1);
    w = zeros(p, 1);
    for n = p+1:N
        xvec = x(n-1:-1:n-p);
        yhat(n) = w.' * xvec;
        e(n) = x(n) - yhat(n);
        w = w + mu * e(n) * xvec;
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
