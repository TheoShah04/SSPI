clear; clc;
dataDir = fullfile('Recordings_27_28_01_20_Shared');
fileName = 'P005.csv';
filePath = fullfile(dataDir, fileName);

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

seg1 = s1:e1;
seg2 = s2:e2;
seg3 = s3:e3;

figure;
plot(ts(seg1), xECG(seg1), 'k');
title('P005 - Segment 1');
xlabel('Time sample');
ylabel('ECG (CH1)');
grid on;

figure;
plot(ts(seg2), xECG(seg2), 'k');
title('P005 - Segment 2');
xlabel('Time sample');
ylabel('ECG (CH1)');
grid on;

figure;
plot(ts(seg3), xECG(seg3), 'k');
title('P005 - Segment 3');
xlabel('Time sample');
ylabel('ECG (CH1)');
grid on;

segments = {xECG(seg1), xECG(seg2), xECG(seg3)};
RRIs = cell(3, 1);
fsRRI = zeros(3, 1);
minLen = fsECG * 10;

hasButter = exist('butter', 'file') == 2;
hasFiltfilt = exist('filtfilt', 'file') == 2;
if ~(hasButter && hasFiltfilt)
    warning('ECG_to_RRI requires the Signal Processing Toolbox. Skipping RRI conversion.');
else
    for k = 1:3
        if numel(segments{k}) < minLen
            warning('Segment %d is shorter than 10 s; skipping RRI.', k);
            continue;
        end
        try
            [RRIs{k}, fsRRI(k)] = ECG_to_RRI(segments{k}, fsECG);
        catch ME
            warning('ECG_to_RRI failed for segment %d: %s', k, ME.message);
        end
    end
end

% ---- Autocorrelation of RRI for three trials ----
maxLag = 50;
for k = 1:3
    if isempty(RRIs{k})
        warning('RRIs{%d} is empty; skipping autocorrelation.', k);
        continue;
    end
    rri = RRIs{k}(:);
    rri = rri(isfinite(rri));
    if numel(rri) < 5
        warning('RRIs{%d} too short; skipping autocorrelation.', k);
        continue;
    end
    rri = detrend(rri, 'constant');
    L = min(maxLag, numel(rri) - 1);
    [acf, lags] = local_autocorr(rri, L);

    figure;
    stem(lags, acf, 'filled');
    title(sprintf('RRI Autocorrelation - Trial %d', k));
    xlabel('Lag (samples)');
    ylabel('Autocorrelation');
    grid on;
end

% ---- AR model order selection for RRI (PACF + AIC/MDL/AICc) ----
maxp = 20;
for trial = 1:3
    if isempty(RRIs{trial})
        warning('RRIs{%d} is empty; skipping AR model selection.', trial);
        continue;
    end
    x = RRIs{trial}(:);
    x = x(isfinite(x));
    if numel(x) < maxp + 5
        warning('Trial %d RRI too short for maxp=%d. Skipping.', trial, maxp);
        continue;
    end

    x = detrend(x, 'constant');
    sx = std(x);
    if sx > 0
        x = x / sx;
    end

    [pacf, conf] = pacf_yw(x, maxp);
    p_pacf = last_significant_lag(pacf, conf);

    figure;
    stem(1:maxp, pacf, 'filled');
    hold on;
    yline(conf, 'r--'); yline(-conf, 'r--');
    hold off;
    xlabel('Lag');
    ylabel('PACF');
    title(sprintf('RRI PACF (Trial %d)', trial));
    grid on;

    [MDL, AIC, AICc] = order_select_yw(x, maxp);
    [~, p_mdl] = min(MDL);
    [~, p_aic] = min(AIC);
    [~, p_aicc] = min(AICc);

    figure;
    plot(1:maxp, MDL, '-o', 'DisplayName', 'MDL');
    hold on;
    plot(1:maxp, AIC, '-s', 'DisplayName', 'AIC');
    plot(1:maxp, AICc, '-d', 'DisplayName', 'AICc');
    hold off;
    xlabel('Model order p');
    ylabel('Criterion value');
    title(sprintf('Order selection (Trial %d)', trial));
    legend('show', 'Location', 'northeast');
    grid on;

    fprintf('Trial %d: PACF suggests p=%d | MDL=%d | AIC=%d | AICc=%d\n', ...
        trial, p_pacf, p_mdl, p_aic, p_aicc);
end

figure;
plot(ts, xECG, 'k');
title('P005');
xlabel('Time sample');
ylabel('ECG (CH1)');
grid on;

% ---- Heart rate from Trial 1 (unconstrained) and PDE plots ----
if ~isempty(RRIs{1})
    rri = RRIs{1}(:);
    rri = rri(~isnan(rri) & rri > 0);
    h = 60 ./ rri;
    h = h(isfinite(h));

    blockSize = 10;
    nBlocks = floor(numel(h) / blockSize);
    hTrim = h(1:blockSize * nBlocks);
    hBlocks = reshape(hTrim, blockSize, []);

    alphas = [1, 0.6];
    for a = alphas
        bh = a * mean(hBlocks, 1);
        figure;
        hold on;
        local_plot_pde(h, 'k', 'h[n]');
        labelStr = ['\hat{h}[n], \alpha=' num2str(a, '%.1f')];
        local_plot_pde(bh, [0 0.6 0], labelStr);
        hold off;
        title(['PDE of Heart Rate (Trial 1) - \alpha=' num2str(a, '%.1f')], "FontSize", 14);
        xlabel('Heart rate (bpm)', 'FontSize', 14);
        ylabel('Probability density', 'FontSize', 14);
        legend('show');
        grid on;
    end
else
    warning('RRIs{1} is empty. Run ECG_to_RRI first to generate Trial 1 RRI.');
end

function idx = local_idx_ge(ts, t, n)
    % First index where ts >= t (clamped to [1, n])
    idx = find(ts >= t, 1, 'first');
    if isempty(idx)
        idx = n;
    end
end

function idx = local_idx_le(ts, t, n)
    % Last index where ts <= t (clamped to [1, n])
    idx = find(ts <= t, 1, 'last');
    if isempty(idx)
        idx = 1;
    end
end

function local_plot_pde(x, colorSpec, labelStr)
    x = x(:);
    x = x(isfinite(x));
    n = numel(x);
    if n < 2
        return;
    end

    sigma = std(x);
    if ~isfinite(sigma) || sigma == 0
        sigma = max(eps, range(x) / 6);
    end
    if exist('ksdensity', 'file') == 2
        bw = 1.06 * sigma * n^(-1/5);
        [f, xi] = ksdensity(x, 'Bandwidth', bw);
        plot(xi, f, 'Color', colorSpec, 'LineWidth', 1.5, 'DisplayName', labelStr);
    else
        binWidth = 3.5 * sigma * n^(-1/3);
        if ~isfinite(binWidth) || binWidth <= 0
            binWidth = max(eps, range(x) / 40);
        end
        histogram(x, 'BinWidth', binWidth, 'Normalization', 'pdf', ...
            'EdgeColor', colorSpec, 'DisplayName', labelStr, ...
            'DisplayStyle', 'stairs', 'LineWidth', 1.2);
    end
end

function [acf, lags] = local_autocorr(x, maxLag)
    x = x(:);
    n = numel(x);
    acf_pos = zeros(maxLag + 1, 1);
    for k = 0:maxLag
        acf_pos(k + 1) = (x(1:n-k)' * x(1+k:n)) / n;
    end
    if acf_pos(1) ~= 0
        acf_pos = acf_pos / acf_pos(1);
    end
    acf = [flipud(acf_pos(2:end)); acf_pos];
    lags = (-maxLag:maxLag).';
end

function [pacf, conf] = pacf_yw(x, maxp)
    x = x(:);
    N = length(x);
    pacf = zeros(maxp, 1);

    for p = 1:maxp
        r0p = acf_unbiased(x, p);
        R = toeplitz(r0p(1:p));   
        rvec = r0p(2:p+1); 
        a = R \ rvec;   
        pacf(p) = a(end);  
    end

    conf = 1.96 / sqrt(N);
end

function p = last_significant_lag(pacf, conf)
    idx = find(abs(pacf) > conf);
    if isempty(idx)
        p = 0;
    else
        p = idx(end);
    end
end

function [MDL, AIC, AICc] = order_select_yw(x, maxp)
    x = x(:);
    N = length(x);
    MDL = zeros(maxp, 1);
    AIC = zeros(maxp, 1);
    AICc = zeros(maxp, 1);

    for p = 1:maxp
        r0p = acf_unbiased(x, p);   
        R = toeplitz(r0p(1:p)); 
        rvec = r0p(2:p+1); 
        a = R \ rvec; 

        e = zeros(N - p, 1);
        for n = p+1:N
            xhat = a' * x(n-1:-1:n-p);
            e(n - p) = x(n) - xhat;
        end

        Ep = sum(e.^2);
        MDL(p) = log(Ep) + (p * log(N)) / N;
        AIC(p) = log(Ep) + (2 * p) / N;
        AICc(p) = AIC(p) + (2 * p * (p + 1)) / (N - p - 1);
    end
end

function r = acf_unbiased(x, maxp)
    x = x(:);
    n = numel(x);
    r = zeros(maxp + 1, 1);
    for k = 0:maxp
        r(k + 1) = (x(1:n-k)' * x(1+k:n)) / (n - k);
    end
end
