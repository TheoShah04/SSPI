clear; clc;

dataDir = fullfile('Recordings_27_28_01_20_Shared');
if ~isfolder(dataDir)
    dataDir = fullfile('task2', 'Recordings_27_28_01_20_Shared');
end
fileName = 'P005.csv';
filePath = fullfile(dataDir, fileName);

if ~isfile(filePath)
    error('Could not find input file: %s', filePath);
end

data = readmatrix(filePath, 'NumHeaderLines', 1);
ts = data(:, 1);
xECG = data(:, 2);
fsECG = 500; % Hz

bounds = [7991400, 15823900];
pad = 100000;
n = numel(xECG);

t1a = bounds(1) - pad;
t1b = bounds(1) + pad;
t2a = bounds(2) - pad;
t2b = bounds(2) + pad;

s1 = 1;
e1 = local_idx_le(ts, t1a, n);
s2 = local_idx_ge(ts, t1b, n);
e2 = local_idx_le(ts, t2a, n);
s3 = local_idx_ge(ts, t2b, n);
e3 = n;

segments = {xECG(s1:e1), xECG(s2:e2), xECG(s3:e3)};
RRIs = cell(3, 1);
fsRRI = zeros(3, 1);

for k = 1:3
    if numel(segments{k}) < fsECG * 10
        warning('Segment %d shorter than 10 s; skipping RRI conversion.', k);
        continue;
    end
    try
        [RRIs{k}, fsRRI(k)] = ECG_to_RRI(segments{k}, fsECG);
    catch ME
        warning('ECG_to_RRI failed for segment %d: %s', k, ME.message);
    end
end

winSecList = [50, 150];
skipBins = 5; % exclude first few low-frequency bins (including DC) in plots

for trial = 1:3
    if isempty(RRIs{trial}) || fsRRI(trial) <= 0
        warning('Trial %d has no valid RRI; skipping PSD.', trial);
        continue;
    end

    x = RRIs{trial}(:);
    x = x(isfinite(x));
    if numel(x) < 16
        warning('Trial %d RRI too short; skipping PSD.', trial);
        continue;
    end
    x = detrend(x, 'constant'); % remove DC before PSD

    [fStd, Pstd] = local_standard_periodogram(x, fsRRI(trial));

    [fAvg50, Pavg50, K50] = local_averaged_periodogram(x, fsRRI(trial), winSecList(1));
    [fAvg150, Pavg150, K150] = local_averaged_periodogram(x, fsRRI(trial), winSecList(2));

    idxStd = (skipBins + 1):numel(fStd);
    figure;
    semilogy(fStd(idxStd), Pstd(idxStd) + eps, 'k', 'LineWidth', 1.2, 'DisplayName', 'Standard periodogram');
    hold on;
    if ~isempty(Pavg50)
        idx50 = (skipBins + 1):numel(fAvg50);
        semilogy(fAvg50(idx50), Pavg50(idx50) + eps, 'b', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('Averaged periodogram (50 s, K=%d)', K50));
    end
    if ~isempty(Pavg150)
        idx150 = (skipBins + 1):numel(fAvg150);
        semilogy(fAvg150(idx150), Pavg150(idx150) + eps, 'r', 'LineWidth', 1.2, ...
            'DisplayName', sprintf('Averaged periodogram (150 s, K=%d)', K150));
    end
    hold off;
    grid on;
    xlabel('Frequency (Hz)', 'FontSize', 14);
    ylabel('PSD', 'FontSize', 14);
    title(sprintf('RRI PSD - Trial %d', trial), 'FontSize', 14);
    legend('show', 'Location', 'northeast');
    set(gca, 'FontSize', 14);
    xlim([0, fsRRI(trial)/2]);
end

function idx = local_idx_ge(ts, t, n)
idx = find(ts >= t, 1, 'first');
if isempty(idx), idx = n; end
end

function idx = local_idx_le(ts, t, n)
idx = find(ts <= t, 1, 'last');
if isempty(idx), idx = 1; end
end

function [f, P1] = local_standard_periodogram(x, fs)
x = x(:);
N = numel(x);
X = fft(x, N);
P2 = (1 / (fs * N)) * abs(X).^2;
nHalf = floor(N/2);
P1 = P2(1:nHalf+1);
if nHalf > 1
    P1(2:end-1) = 2 * P1(2:end-1);
end
f = (0:nHalf)' * (fs / N);
end

function [f, Pavg, K] = local_averaged_periodogram(x, fs, winSec)
x = x(:);
N = numel(x);
L = round(winSec * fs);% window length in samples
if L < 8 || N < L
    f = [];
    Pavg = [];
    K = 0;
    return;
end

K = floor(N / L);% non-overlapping segments
if K < 1
    f = [];
    Pavg = [];
    return;
end

w = hann(L, 'periodic');
U = mean(w.^2);% window power normalisation
nHalf = floor(L/2);
Psum = zeros(nHalf+1, 1);

for k = 1:K
    idx = (k-1)*L + (1:L);
    xk = x(idx) .* w;
    Xk = fft(xk, L);
    P2 = (1 / (fs * L * U)) * abs(Xk).^2;
    P1 = P2(1:nHalf+1);
    if nHalf > 1
        P1(2:end-1) = 2 * P1(2:end-1);
    end
    Psum = Psum + P1;
end

Pavg = Psum / K;
f = (0:nHalf)' * (fs / L);
end
