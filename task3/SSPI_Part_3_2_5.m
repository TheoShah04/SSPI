clear; clc;

orders = [2, 10, 30];   % low, medium, high model order
Nfft = 1024;

% Load sunspot series (year, sunspot number)
S = load('sunspot.dat');
if isstruct(S)
    fn = fieldnames(S);
    D = S.(fn{1});
else
    D = S;
end

if size(D, 2) >= 2
    year = D(:, 1);
    sunspot = D(:, 2);
else
    year = (1:numel(D)).';
    sunspot = D(:);
end

sunspot = sunspot(:);
year = year(:);

% Original and mean-centered versions
series = {sunspot, sunspot - mean(sunspot)};
seriesLabel = {'Original series', 'Mean-centered series'};

for s = 1:numel(series)
    x = series{s};
    N = numel(x);

    % One-sided periodogram using pgm.m
    Pp = pgm(x);
    f = (0:N-1) / N;
    idx = 1:(floor(N/2) + 1);

    figure;
    t = tiledlayout(1, numel(orders), 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t, sprintf('Sunspot PSD: %s', seriesLabel{s}), 'FontSize', 16);

    for k = 1:numel(orders)
        p = orders(k);
        [aHat, sigma2Hat] = aryule(x, p);
        [h, w] = freqz(sqrt(sigma2Hat), aHat, Nfft);
        Pm = abs(h).^2;

        nexttile;
        semilogy(f(idx), Pp(idx), 'Color', [0.6 0.6 0.6], 'DisplayName', 'Periodogram');
        hold on;
        semilogy(w/(2*pi), Pm, 'b', 'LineWidth', 1.3, ...
            'DisplayName', sprintf('Model-based AR(%d)', p));
        hold off;
        grid on;
        xlim([0 0.5]);
        xlabel('Normalized frequency', 'FontSize', 14);
        ylabel('PSD', 'FontSize', 14);
        title(sprintf('AR(%d)', p), 'FontSize', 14);
        set(gca, 'FontSize', 14);
        if k == 1
            legend('show', 'Location', 'northeast', 'FontSize', 12);
        end
    end
end

% Time-series plots for reference
figure;
subplot(2,1,1);
plot(year, sunspot, 'k');
grid on;
xlabel('Year', 'FontSize', 14);
ylabel('Sunspot number', 'FontSize', 14);
title('Sunspot series (original)', 'FontSize', 14);
set(gca, 'FontSize', 14);

subplot(2,1,2);
plot(year, sunspot - mean(sunspot), 'b');
grid on;
xlabel('Year', 'FontSize', 14);
ylabel('Sunspot number', 'FontSize', 14);
title('Sunspot series (mean-centered)', 'FontSize', 14);
set(gca, 'FontSize', 14);

fprintf('\nInterpretation guide:\n');
fprintf('- Original series includes strong DC/very-low-frequency energy from non-zero mean.\n');
fprintf('- Mean-centering removes DC bias and reveals oscillatory content more clearly.\n');
fprintf('- Under-modelling (small AR order) oversmooths PSD and misses narrow/secondary peaks.\n');
fprintf('- Over-modelling (large AR order) can create spurious sharp peaks and fit noise.\n');
fprintf('- Moderate order gives a smoother PSD that captures dominant spectral structure.\n\n');
