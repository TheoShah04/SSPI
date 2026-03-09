clear; clc; close all;

N = 1000;
sigma = 0.1;
b = [1, 2, 3, 2, 1];
a = 1;
Nw = numel(b) - 1;              % Adaptive filter order (5 coefficients total)
M = Nw + 1;

rng(42);

x = randn(N, 1);
y_raw = filter(b, a, x);

% Normalize y to unit variance
scale = 1 / std(y_raw, 1);
y = scale * y_raw;

eta = sigma * randn(N, 1);
z = y + eta;

% Wiener solution estimated from sample correlations
lagsRxx = -Nw:Nw;
rxx = corr_est(x, x, lagsRxx);
idx0 = find(lagsRxx == 0, 1);
rxx_pos = rxx(idx0:end);
rxx_neg = rxx(idx0:-1:1);
Rxx = toeplitz(rxx_pos, rxx_neg);

lagsRzx = -Nw:0;
rzx = corr_est(z, x, lagsRzx);
pzx = rzx(end:-1:1).';
wWiener = Rxx \ pzx;

lambdaMax = max(real(eig(Rxx)));
muMaxTheory = 2 / lambdaMax;
Px = mean(x .^ 2);
muMaxPractical = 2 / ((Nw + 1) * Px);

% Requested case: mu = 0.01
mu0 = 0.01;
[~, e0, W0] = lms(x, z, mu0, Nw);

fprintf('LMS study (N=%d, sigma=%.3f, Nw=%d)\n', N, sigma, Nw);
fprintf('Estimated mean-convergence bound: 0 < mu < %.4f\n', muMaxTheory);
fprintf('Practical mean-square bound (white-input approx): mu < %.4f\n', muMaxPractical);
fprintf('Requested run: mu = %.4f\n', mu0);
fprintf('Wiener coefficients wW^T = ');
fprintf('%.4f ', wWiener);
fprintf('\n');
fprintf('Final LMS coefficients w_end^T = ');
fprintf('%.4f ', W0(:, end));
fprintf('\n\n');

% Plot all 5 coefficient trajectories on one graph for mu=0.01
cmap = lines(M);
figure('Name', 'Coefficient Evolution (\mu = 0.01)');
hold on;
for k = 1:M
    plot(W0(k, :), 'Color', cmap(k, :), 'LineWidth', 1.1, ...
        'DisplayName', sprintf('w_%d', k - 1));
    yline(wWiener(k), '--', 'Color', cmap(k, :), 'HandleVisibility', 'off');
end
hold off;
grid on;
xlabel('n');
ylabel('Coefficient value');
title('LMS coefficient evolution (\mu = 0.01)');
legend('Location', 'best');

% Show squared estimate error for mu=0.01
figure('Name', 'Squared Estimate Error (\mu = 0.01)');
plot(e0 .^ 2, 'Color', [0.15 0.15 0.15], 'LineWidth', 0.8, 'DisplayName', 'e^2[n]');
hold on;
plot(movmean(e0 .^ 2, 50), 'r', 'LineWidth', 1.2, ...
    'DisplayName', '50-sample moving average');
hold off;
grid on;
xlabel('n');
ylabel('Squared error');
title('Squared estimate error (\mu = 0.01)');
legend('Location', 'best');

% Sweep mu in [0.002, 0.5]
muVals = [0.002, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50];
numMu = numel(muVals);

W_all = cell(numMu, 1);
e2smooth_all = cell(numMu, 1);
steadyMSE = zeros(numMu, 1);
finalRelErr = zeros(numMu, 1);
convSample = nan(numMu, 1);
isDiverged = false(numMu, 1);

tol = 0.05;         % 5% relative distance to Wiener solution
holdSamples = 50;   % must remain below threshold for this many samples

for i = 1:numMu
    mu = muVals(i);
    [~, ei, Wi] = lms(x, z, mu, Nw);

    e2 = ei .^ 2;
    relDist = sqrt(sum((Wi - wWiener) .^ 2, 1)) / max(norm(wWiener), eps);
    wThreshold = 25 * max(1, max(abs(wWiener)));
    isDiverged(i) = any(~isfinite(Wi(:))) || any(abs(Wi(:)) > wThreshold) || ...
        any(~isfinite(e2)) || any(e2 > 1e4);

    W_all{i} = Wi;
    e2smooth_all{i} = movmean(e2, 50);
    steadyMSE(i) = mean(e2(round(0.8 * N):end));
    finalRelErr(i) = relDist(end);
    convSample(i) = find_conv_index(relDist, tol, holdSamples);
end

% Plot all 5 coefficients together for each mu (one tile per mu)
figure('Name', 'Coefficient Evolution for \mu Sweep');
tl = tiledlayout(4, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:numMu
    nexttile;
    hold on;
    Wi = W_all{i};
    for k = 1:M
        plot(Wi(k, :), 'Color', cmap(k, :), 'LineWidth', 0.9);
        yline(wWiener(k), '--', 'Color', cmap(k, :), 'HandleVisibility', 'off');
    end
    hold off;
    grid on;
    xlabel('n', 'FontSize', 16);
    ylabel('w', 'FontSize', 16);
    if isDiverged(i)
        title(sprintf('\\mu = %.3f (unstable/outlier)', muVals(i)), 'FontSize', 14);
    else
        title(sprintf('\\mu = %.3f', muVals(i)), 'FontSize', 14);
    end
    set(gca, 'FontSize', 14);
end
sgtitle(tl, 'Coefficient Evolution for Different \mu', 'FontSize', 18);

% Smoothed squared error for different mu values
figure('Name', 'Smoothed Squared Error for \mu Sweep');
hold on;
for i = 1:numMu
    semilogy(max(e2smooth_all{i}, eps), 'LineWidth', 1.1, ...
        'DisplayName', sprintf('\\mu = %.3f', muVals(i)));
end
hold off;
grid on;
xlabel('n');
ylabel('Smoothed e^2[n] (log scale)');
title('Learning curves for different adaptation gains (log y-axis)');
legend('Location', 'northeast');

% Optional zoom: stable-only learning curves for readability
if any(~isDiverged)
    figure('Name', 'Smoothed Squared Error (\mu stable subset)');
    hold on;
    for i = 1:numMu
        if ~isDiverged(i)
            plot(e2smooth_all{i}, 'LineWidth', 1.1, ...
                'DisplayName', sprintf('\\mu = %.3f', muVals(i)));
        end
    end
    hold off;
    grid on;
    xlabel('n');
    ylabel('Smoothed e^2[n]');
    title('Learning curves excluding unstable/outlier \mu values');
    legend('Location', 'northeast');
end

% Summary metrics vs mu
figure('Name', 'Convergence Metrics vs \mu');
tiledlayout(2, 1);
nexttile;
muPlotMask = muVals < 0.5;
loglog(muVals(muPlotMask), max(steadyMSE(muPlotMask), eps), 'o-', 'LineWidth', 1.2);
grid on;
xlabel('\mu');
ylabel('Steady-state mean(e^2)');
title('Steady-state error vs adaptation gain (log-log, \mu < 0.5)');

nexttile;
semilogx(muVals, convSample, 's-', 'LineWidth', 1.2);
grid on;
xlabel('\mu');
ylabel('Convergence sample index');
title(sprintf('Convergence to Wiener (< %.0f%% for %d samples)', 100 * tol, holdSamples));

function [yhat, e, W] = lms(x, z, mu, Nw)
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

function idx = find_conv_index(relDist, tol, holdSamples)
    idx = NaN;
    N = numel(relDist);
    if N < holdSamples
        return;
    end

    for n = 1:(N - holdSamples + 1)
        if all(relDist(n:n + holdSamples - 1) < tol)
            idx = n;
            return;
        end
    end
end
