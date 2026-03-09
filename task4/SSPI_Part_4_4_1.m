clear; clc; close all;

% AR(2) synthesis model in MATLAB notation:
% x[n] + 0.9 x[n-1] + 0.2 x[n-2] = eta[n]
% => x[n] = -0.9 x[n-1] - 0.2 x[n-2] + eta[n]

N = 4000;
a = [1, 0.9, 0.2];
muVals = [0.005, 0.01, 0.03, 0.05];   % Four adaptation gains (includes 0.01)

rng(32);
eta = randn(N, 1);
x = filter(1, a, eta);

% Expected predictor coefficients from AR sign convention
wExpected = -a(2:end).';   % [a1* a2*] = [-0.9 -0.2]
M = numel(wExpected);

numMu = numel(muVals);
W_all = cell(numMu, 1);
e_all = cell(numMu, 1);
wFinal = zeros(M, numMu);
relErr = zeros(numMu, 1);
ssMSE = zeros(numMu, 1);
convIdxA1 = nan(numMu, 1);
convIdxA2 = nan(numMu, 1);
jitterA1 = zeros(numMu, 1);
jitterA2 = zeros(numMu, 1);

ssWin = 800;              % steady-state window length
tolPct = 0.10;            % 10% of target value
holdSamples = 40;         % evaluation window length for settling
requiredFrac = 0.80;      % at least 80% of samples in window must be in tolerance
tailWin = 800;            % jitter window

for i = 1:numMu
    mu = muVals(i);
    [W, e] = lms_ar2_predictor(x, mu);

    W_all{i} = W;
    e_all{i} = e;
    wFinal(:, i) = W(:, end);
    relErr(i) = norm(wFinal(:, i) - wExpected) / norm(wExpected);
    ssMSE(i) = mean(e(end - ssWin + 1:end) .^ 2);

    tolA1 = tolPct * abs(wExpected(1));
    tolA2 = tolPct * abs(wExpected(2));
    convIdxA1(i) = settle_index(W(1, :), wExpected(1), tolA1, holdSamples, requiredFrac);
    convIdxA2(i) = settle_index(W(2, :), wExpected(2), tolA2, holdSamples, requiredFrac);

    a1TailErr = W(1, end - tailWin + 1:end) - wExpected(1);
    a2TailErr = W(2, end - tailWin + 1:end) - wExpected(2);
    jitterA1(i) = std(a1TailErr, 1);
    jitterA2(i) = std(a2TailErr, 1);
end

fprintf('AR model (MATLAB): a = [%.2f %.2f %.2f]\n', a);
fprintf('Expected LMS predictor coefficients: [a1 a2] = [%.4f %.4f]\n', ...
    wExpected(1), wExpected(2));
fprintf(['Convergence criterion per coefficient: |a_k[n]-a_k*| <= %.1f%%%%|a_k*| ', ...
    'for at least %.0f%%%% of a %d-sample window\n\n'], 100 * tolPct, 100 * requiredFrac, holdSamples);
fprintf('%-8s %-10s %-10s %-10s %-12s %-10s %-10s %-10s %-10s\n', ...
    'mu', 'a1(final)', 'a2(final)', 'RelErr', 'SS-MSE', 'nConv(a1)', 'nConv(a2)', 'jit(a1)', 'jit(a2)');

for i = 1:numMu
    if isnan(convIdxA1(i)), c1 = 'n/a'; else, c1 = sprintf('%d', round(convIdxA1(i))); end
    if isnan(convIdxA2(i)), c2 = 'n/a'; else, c2 = sprintf('%d', round(convIdxA2(i))); end
    fprintf('%-8.3f %-10.4f %-10.4f %-10.4f %-12.6f %-10s %-10s %-10.4f %-10.4f\n', ...
        muVals(i), wFinal(1, i), wFinal(2, i), relErr(i), ssMSE(i), c1, c2, jitterA1(i), jitterA2(i));
end

fprintf('\nInterpretation guide:\n');
fprintf('- Smaller mu -> slower convergence (larger nConv), but lower jitter in steady state.\n');
fprintf('- Larger mu -> faster convergence (smaller nConv), but larger coefficient fluctuations and SS-MSE.\n');
fprintf('- For this AR(2) setup, a1 often settles earlier than a2 because lag-1 structure is stronger.\n');

% Plot 1: Coefficient evolution for each mu
figure('Name', 'AR(2) LMS Coefficient Evolution for Different \mu');
tiledlayout(2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for i = 1:numMu
    nexttile;
    W = W_all{i};
    plot(W(1, :), 'LineWidth', 1.2, 'DisplayName', 'a_1[n]');
    hold on;
    plot(W(2, :), 'LineWidth', 1.2, 'DisplayName', 'a_2[n]');
    yline(wExpected(1), '--', 'HandleVisibility', 'off');
    yline(wExpected(2), '--', 'HandleVisibility', 'off');
    hold off;
    grid on;
    xlabel('n', 'FontSize', 14);
    ylabel('Coefficient value', 'FontSize', 14);
    title(sprintf('\\mu = %.3f', muVals(i)), 'FontSize', 14);
    if i == 1
        legend('Location', 'best');
    end
    set(gca, 'FontSize', 14);
end

% Plot 2: Squared-error learning curves for each mu
figure('Name', 'AR(2) LMS Squared Error for Different \mu');
hold on;
for i = 1:numMu
    e = e_all{i};
    semilogy(movmean(e .^ 2, 100), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('\\mu = %.3f', muVals(i)));
end
hold off;
grid on;
xlabel('n', 'FontSize', 14);
ylabel('Smoothed e^2[n] (log)', 'FontSize', 14);
title('Prediction Error vs n for Different \mu', 'FontSize', 14);
legend('Location', 'best');
set(gca, 'FontSize', 14);

% Plot 3: Convergence in (a1, a2) parameter space
figure('Name', 'AR(2) Parameter-Space Convergence for Different \mu');
hold on;
for i = 1:numMu
    W = W_all{i};
    plot(W(1, 3:end), W(2, 3:end), 'LineWidth', 1.2, ...
        'DisplayName', sprintf('\\mu = %.3f', muVals(i)));
    plot(W(1, end), W(2, end), 'o', 'HandleVisibility', 'off');
end
plot(wExpected(1), wExpected(2), 'rp', 'MarkerSize', 12, 'MarkerFaceColor', 'r', ...
    'DisplayName', 'expected');
hold off;
grid on;
xlabel('a_1', 'FontSize', 14);
ylabel('a_2', 'FontSize', 14);
title('Coefficient Trajectories in Parameter Space', 'FontSize', 14);
legend('Location', 'best');
axis equal;
set(gca, 'FontSize', 14);

function [W, e] = lms_ar2_predictor(x, mu)
    N = numel(x);
    W = zeros(2, N);
    e = zeros(N, 1);
    w = zeros(2, 1);

    for n = 3:N
        xvec = [x(n - 1); x(n - 2)];
        yhat = w.' * xvec;
        e(n) = x(n) - yhat;
        w = w + mu * e(n) * xvec;
        W(:, n) = w;
    end
end

function idx = settle_index(track, target, tolAbs, holdSamples, requiredFrac)
    idx = NaN;
    if numel(track) < holdSamples
        return;
    end

    inTol = abs(track - target) <= tolAbs;
    for n = 1:(numel(track) - holdSamples + 1)
        if mean(inTol(n:n + holdSamples - 1)) >= requiredFrac
            idx = n;
            return;
        end
    end
end
