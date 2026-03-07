clear; clc;

h = 0.2 * ones(1, 5);
groupDelay = 2; % (L-1)/2 for L=5

rng(42); % reproducible
Ns = [128, 256, 512];

for i = 1:numel(Ns)
    N = Ns(i);
    x = randn(N, 1);

    P = pgm(x);
    f = (0:N-1) / N;

    Psm = filter(h, 1, P);
    Psm = circshift(Psm, -groupDelay); % align for zero-phase display

    figure;
    plot(f, P, 'Color', [0.6 0.6 0.6], 'DisplayName', 'Periodogram');
    hold on;
    plot(f, Psm, 'b', 'LineWidth', 1.2, 'DisplayName', 'Smoothed (MA-5)');
    hold off;
    title(sprintf('Periodogram Smoothing (N=%d)', N), 'FontSize', 16);
    xlabel('Normalised frequency', 'FontSize', 14);
    ylabel('P_b_x(f)', 'FontSize', 14);
    legend('show', 'Location', 'northeast');
    grid on;
    set(gca, 'FontSize', 14);

    fprintf('N=%d: mean(Psm)=%.4f, var(Psm)=%.4f\n', N, mean(Psm), var(Psm));
end
