clear; clc; close all;

N = 1000;
sigma = 0.1;
b = [1, 2, 3, 2, 1];
a = 1;
Nw = numel(b) - 1;
M = Nw + 1;

rng(42);
x = randn(N, 1);
y_raw = filter(b, a, x);
scale = 1 / std(y_raw, 1);
y = scale * y_raw;
wTrue = (scale * b).';
eta = sigma * randn(N, 1);
z = y + eta;

% Wiener benchmark (sample estimate)
lagsRxx = -Nw:Nw;
rxx = corr_est(x, x, lagsRxx);
idx0 = find(lagsRxx == 0, 1);
Rxx = toeplitz(rxx(idx0:end), rxx(idx0:-1:1));

lagsRzx = -Nw:0;
rzx = corr_est(z, x, lagsRzx);
pzx = rzx(end:-1:1).';
wWiener = Rxx \ pzx;

lambdaMax = max(real(eig(Rxx)));
muTheoryMax = 2 / lambdaMax;
Px = mean(x .^ 2);
muPracticalMax = 2 / ((Nw + 1) * Px);

% Fixed-mu baseline for comparison
muFixed = 0.01;
[~, eFixed, WFixed] = lms_fixed(x, z, muFixed, Nw);

% Gear-shifting settings
cfg.muMin = 0.002;
cfg.muMax = min(0.25, 0.8 * muPracticalMax);   
cfg.muMax = max(cfg.muMax, 2 * cfg.muMin);    
cfg.alphaErr = 0.98;                       
cfg.betaMu = 0.90;                        
cfg.errLow = sigma^2;                
cfg.errHigh = 10 * sigma^2;  

[~, eGear, WGear, muHist, errPowHist] = lms_gear(x, z, Nw, cfg);

fprintf('Gear-shifting LMS (N=%d, M=%d)\n', N, M);
fprintf('Mean-convergence bound:      0 < mu < %.4f\n', muTheoryMax);
fprintf('Practical LMS bound (approx):0 < mu < %.4f\n', muPracticalMax);
fprintf('Chosen mu range:             [%.4f, %.4f]\n', cfg.muMin, cfg.muMax);

% Quantitative comparison
win = 200;
mseFixed_ss = mean(eFixed(end-win+1:end) .^ 2);
mseGear_ss = mean(eGear(end-win+1:end) .^ 2);
fprintf('Steady-state mean(e^2), fixed mu=%.3f: %.6f\n', muFixed, mseFixed_ss);
fprintf('Steady-state mean(e^2), gear shifting: %.6f\n', mseGear_ss);

% Overshoot/rise-time comparison data (fixed mu=0.01 vs gear shifting only)
numMethods = 2;
methodNames = {'fixed_0.010'; 'gear_shift'};
WMethods = {WFixed; WGear};
ssMSE = [mseFixed_ss; mseGear_ss];

overshootLimitPct = 20;
overshootPct = zeros(numMethods, M);
riseTimes = nan(numMethods, M);
passOvershoot = false(numMethods, 1);
riseWorst = nan(numMethods, 1);
riseMean = nan(numMethods, 1);

for i = 1:numMethods
    [overshootPct(i, :), riseTimes(i, :)] = coeff_overshoot_risetime(WMethods{i}, wTrue);
    passOvershoot(i) = all(overshootPct(i, :) <= overshootLimitPct);

    validRise = isfinite(riseTimes(i, :));
    if all(validRise)
        riseWorst(i) = max(riseTimes(i, :));
        riseMean(i) = mean(riseTimes(i, :));
    end
end

fprintf('True coefficients (scaled): ');
fprintf('%.4f ', wTrue);
fprintf('\n');
fprintf('%-14s %-8s %-11s %-11s %-11s\n', 'Method', 'Pass20%', 'WorstRise', 'MeanRise', 'SS-MSE');
for i = 1:numMethods
    if passOvershoot(i), passTxt = 'yes'; else, passTxt = 'no'; end
    if isnan(riseWorst(i)), worstTxt = 'n/a'; else, worstTxt = sprintf('%d', round(riseWorst(i))); end
    if isnan(riseMean(i)), meanTxt = 'n/a'; else, meanTxt = sprintf('%.1f', riseMean(i)); end
    fprintf('%-14s %-8s %-11s %-11s %-11.6f\n', methodNames{i}, passTxt, worstTxt, meanTxt, ssMSE(i));
end

for i = 1:numMethods
    fprintf('\n%s:\n', methodNames{i});
    fprintf('  Overshoot %% per coeff: ');
    fprintf('%.2f ', overshootPct(i, :));
    fprintf('\n');
    fprintf('  Rise time per coeff:   ');
    for k = 1:M
        if isnan(riseTimes(i, k))
            fprintf('n/a ');
        else
            fprintf('%d ', round(riseTimes(i, k)));
        end
    end
    fprintf('\n');
end

validBest = passOvershoot & isfinite(riseWorst);
if any(validBest)
    idxValid = find(validBest);
    [~, localBest] = min(riseWorst(validBest));
    bestIdx = idxValid(localBest);
    fprintf('\nBest method under 20%% overshoot constraint: %s\n', methodNames{bestIdx});
    fprintf('Worst-coefficient rise time = %d samples\n', round(riseWorst(bestIdx)));
else
    fprintf('\nNo method satisfied the 20%% overshoot constraint for all coefficients.\n');
end

% Plot 1: Weight trajectories (gear-shifting) + Wiener references
cmap = lines(M);
figure('Name', 'Gear-Shifting LMS Weights');
hold on;
for k = 1:M
    plot(WGear(k, :), 'Color', cmap(k, :), 'LineWidth', 1.0, ...
        'DisplayName', sprintf('w_%d', k - 1));
    yline(wWiener(k), '--', 'Color', cmap(k, :), 'HandleVisibility', 'off');
end
hold off;
grid on;
xlabel('n', 'FontSize', 13);
ylabel('Weight value', 'FontSize', 13);
title('Time-Varying-\mu LMS: Weight Evolution', 'FontSize', 15);
legend('Location', 'best');
set(gca, 'FontSize', 12);

% Plot 2: Error behaviour (fixed mu vs gear shifting)
figure('Name', 'Error Comparison');
tiledlayout(2, 1);

nexttile;
semilogy(movmean(eFixed .^ 2, 50), 'k--', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Fixed \\mu = %.3f', muFixed));
hold on;
semilogy(movmean(eGear .^ 2, 50), 'b', 'LineWidth', 1.3, ...
    'DisplayName', 'Gear-shifting \mu[n]');
hold off;
grid on;
xlabel('n');
ylabel('Smoothed e^2[n] (log)');
title('Squared Error (moving average)');
legend('Location', 'best');

nexttile;
plot(eFixed, 'Color', [0.6 0.6 0.6], 'LineWidth', 0.7, ...
    'DisplayName', sprintf('Fixed \\mu = %.3f', muFixed));
hold on;
plot(eGear, 'b', 'LineWidth', 0.9, 'DisplayName', 'Gear-shifting \mu[n]');
hold off;
grid on;
xlabel('n');
ylabel('e[n]');
title('Instantaneous Error');
legend('Location', 'best');

% Plot 3: Adaptation-gain control signals
figure('Name', 'Adaptation Gain Control');
tiledlayout(2, 1);

nexttile;
plot(muHist, 'LineWidth', 1.2);
hold on;
yline(cfg.muMin, '--', 'HandleVisibility', 'off');
yline(cfg.muMax, '--', 'HandleVisibility', 'off');
hold off;
grid on;
xlabel('n');
ylabel('\mu[n]');
title('Time-varying adaptation gain');

nexttile;
plot(errPowHist, 'LineWidth', 1.2);
hold on;
yline(cfg.errLow, '--', 'DisplayName', 'errLow');
yline(cfg.errHigh, '--', 'DisplayName', 'errHigh');
hold off;
grid on;
xlabel('n');
ylabel('Smoothed e^2');
title('Error-power criterion used to manipulate \mu[n]');
legend('Location', 'best');

function [yhat, e, W] = lms_fixed(x, z, mu, Nw)
    N = numel(x);
    M = Nw + 1;
    yhat = zeros(N, 1);
    e = zeros(N, 1);
    W = zeros(M, N);
    w = zeros(M, 1);

    for n = 1:N
        xvec = zeros(M, 1);
        for m = 0:Nw
            if n - m >= 1
                xvec(m + 1) = x(n - m);
            end
        end

        yhat(n) = w.' * xvec;
        e(n) = z(n) - yhat(n);
        w = w + mu * e(n) * xvec;
        W(:, n) = w;
    end
end

function [yhat, e, W, muHist, errPowHist] = lms_gear(x, z, Nw, cfg)
    N = numel(x);
    M = Nw + 1;
    yhat = zeros(N, 1);
    e = zeros(N, 1);
    W = zeros(M, N);
    muHist = zeros(N, 1);
    errPowHist = zeros(N, 1);

    w = zeros(M, 1);
    muNow = cfg.muMax;
    errPowNow = cfg.errHigh;

    for n = 1:N
        xvec = zeros(M, 1);
        for m = 0:Nw
            if n - m >= 1
                xvec(m + 1) = x(n - m);
            end
        end

        yhat(n) = w.' * xvec;
        e(n) = z(n) - yhat(n);

        % Subroutine call:
        % compute mu(n) for the (n+1)-coefficient update from e(n).
        [muNow, errPowNow] = update_mu_from_error(e(n), muNow, errPowNow, cfg);

        w = w + muNow * e(n) * xvec;
        W(:, n) = w;
        muHist(n) = muNow;
        errPowHist(n) = errPowNow;
    end
end

function [muNext, errPowNext] = update_mu_from_error(eNow, muPrev, errPowPrev, cfg)
    % Smoothed error power estimate from current error sample
    errPowNext = cfg.alphaErr * errPowPrev + (1 - cfg.alphaErr) * (eNow ^ 2);

    % Piecewise-linear gear-shifting rule
    if errPowNext <= cfg.errLow
        muTarget = cfg.muMin;
    elseif errPowNext >= cfg.errHigh
        muTarget = cfg.muMax;
    else
        ratio = (errPowNext - cfg.errLow) / (cfg.errHigh - cfg.errLow);
        muTarget = cfg.muMin + ratio * (cfg.muMax - cfg.muMin);
    end

    % Smooth and clip mu to avoid abrupt/unstable changes
    muNext = cfg.betaMu * muPrev + (1 - cfg.betaMu) * muTarget;
    muNext = min(max(muNext, cfg.muMin), cfg.muMax);
end

function r = corr_est(x, y, kvec)
    N = numel(x);
    r = zeros(size(kvec));
    for ii = 1:numel(kvec)
        k = kvec(ii);
        n0 = max(1, 1 - k);
        n1 = min(N, N - k);
        idx = n0:n1;
        r(ii) = mean(x(idx) .* y(idx + k));
    end
end

function [overshootPct, riseSamples] = coeff_overshoot_risetime(W, wTrue)
    M = numel(wTrue);
    overshootPct = zeros(1, M);
    riseSamples = nan(1, M);

    for k = 1:M
        wk = W(k, :);
        target = wTrue(k);
        if abs(target) < eps
            overshootPct(k) = NaN;
            riseSamples(k) = NaN;
            continue;
        end

        if target >= 0
            peak = max(wk);
            overshootPct(k) = max(0, (peak - target) / abs(target) * 100);
            n10 = find(wk >= 0.1 * target, 1, 'first');
            n90 = find(wk >= 0.9 * target, 1, 'first');
        else
            trough = min(wk);
            overshootPct(k) = max(0, (target - trough) / abs(target) * 100);
            n10 = find(wk <= 0.1 * target, 1, 'first');
            n90 = find(wk <= 0.9 * target, 1, 'first');
        end

        if ~isempty(n10) && ~isempty(n90) && n90 >= n10
            riseSamples(k) = n90 - n10;
        end
    end
end
