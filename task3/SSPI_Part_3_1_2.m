clear; clc;
rng(42); 
N = 1024;
L = 128;
K = 8;

x = randn(N, 1);
X = reshape(x, L, K); % each column is a segment

% Periodograms for each segment
P = zeros(L, K);
for k = 1:K
    P(:, k) = pgm(X(:, k));
end
f = (0:L-1) / L;

figure;
for k = 1:K
    subplot(2, 4, k);
    plot(f, P(:, k), 'k');
    title(sprintf('Segment %d', k));
    xlabel('Norm. freq');
    ylabel('P_b_x(f)');
    grid on;
end

% Variation across segments
Pmean = mean(P, 2);
Pvar = var(P, 0, 2);

figure;
plot(f, Pmean, 'b', 'LineWidth', 1.2, 'DisplayName', 'Mean');
hold on;
plot(f, Pmean + sqrt(Pvar), 'r--', 'DisplayName', 'Mean \pm Std');
plot(f, Pmean - sqrt(Pvar), 'r--', 'HandleVisibility', 'off');
hold off;
title('Mean Periodogram Across 8 Segments');
xlabel('Normalized frequency');
ylabel('P_b_x(f)');
legend('show', 'Location', 'northeast');
grid on;

Pavg = Pmean;
figure;
plot(f, Pavg, 'b', 'LineWidth', 1.4);
title('Averaged Periodogram (8 Segments)');
xlabel('Normalized frequency');
ylabel('P_b_x(f)');
grid on;

% Segment-wise summary
seg_mean = mean(P, 1);
seg_var = var(P, 0, 1);
disp(table((1:K)', seg_mean', seg_var', ...
    'VariableNames', {'Segment', 'MeanP', 'VarP'}));

fprintf('Mean of segment means: %.4f\n', mean(seg_mean));
fprintf('Mean of segment variances: %.4f\n', mean(seg_var));

figure;
plot(1:K, seg_mean, '-o', 'LineWidth', 1.2, 'DisplayName', 'Mean');
hold on;
plot(1:K, seg_var, '-s', 'LineWidth', 1.2, 'DisplayName', 'Variance');
hold off;
title('Segment-wise Mean and Variance of Periodogram');
xlabel('Segment');
ylabel('Value');
legend('show', 'Location', 'northeast');
grid on;
