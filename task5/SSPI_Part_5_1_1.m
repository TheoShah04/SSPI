clear; clc; close all;

N = 10;
n = (0:N-1).';

% Frequency sweep for true f0
f0Vals = linspace(0.02, 0.48, 9);

fGrid = linspace(0, 0.5, 401);

fHat = zeros(size(f0Vals));
hPlot = gobjects(numel(f0Vals), 1);

figure('Name', 'Periodograms and MLE, N=10');
hold on;

for k = 1:numel(f0Vals)
    f0 = f0Vals(k);
    x = cos(2 * pi * f0 * n);

    % Periodogram on fine grid
    Xf = exp(-1j * 2 * pi * (n * fGrid))' * x; 
    Pxx = (1 / N) * abs(Xf) .^ 2;

    % MLE = argmax of periodogram
    [~, idx] = max(Pxx);
    fHat(k) = fGrid(idx);

    hPlot(k) = plot(fGrid, Pxx, 'LineWidth', 1.0);
    plot(fHat(k), Pxx(idx), 'o', 'MarkerSize', 6, 'HandleVisibility', 'off');
end
hold off;
grid on;
xlim([0 0.5]);
xlabel('f', 'FontSize', 14);
ylabel('P_{XX}(f)', 'FontSize', 14);
title('Periodograms and MLE Peaks (N=10)', 'FontSize', 14);
legLabels = arrayfun(@(f0, fh) sprintf('$f_0=%.3f,\\ \\hat{f}_0=%.3f$', f0, fh), ...
    f0Vals, fHat, 'UniformOutput', false);
lgd = legend(hPlot, legLabels, 'Location', 'best');
set(lgd, 'Interpreter', 'latex');

f0Dense = linspace(0.02, 0.48, 41);
fHatDense = zeros(size(f0Dense));
for k = 1:numel(f0Dense)
    f0 = f0Dense(k);
    x = cos(2 * pi * f0 * n);
    Xf = exp(-1j * 2 * pi * (n * fGrid))' * x;
    Pxx = (1 / N) * abs(Xf) .^ 2;
    [~, idx] = max(Pxx);
    fHatDense(k) = fGrid(idx);
end

figure('Name', 'MLE Error vs f0 (N=10)');
plot(f0Dense, fHatDense - f0Dense, 'o-', 'LineWidth', 1.2);
grid on;
xlabel('f_0', 'FontSize', 14);
ylabel('$\hat{f}_0 - f_0$', 'FontSize', 14, 'Interpreter', 'latex');
title('MLE Frequency Error vs True f_0', 'FontSize', 14);

fprintf('N = %d\n', N);
fprintf('True f0 sweep: [%.3f, %.3f]\n', f0Vals(1), f0Vals(end));
fprintf('Average abs error: %.6f\n', mean(abs(fHat - f0Vals)));
fprintf('Max abs error: %.6f\n', max(abs(fHat - f0Vals)));
