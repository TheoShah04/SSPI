clear; clc;

N_vals = 1:50:1001;
sigma2_vals = 1:50:1001;

[Ngrid, Sgrid] = meshgrid(N_vals, sigma2_vals);

% CRLB for sigma^2: var(sigma^2_hat) >= 2*sigma^4 / N
CRLB_sigma2 = 2 .* (Sgrid .^ 2) ./ Ngrid;

% CRLB for a1: var(a1_hat) >= (1 - a1^2) / N
a1 = -0.005581;
% a1 = 0.979;
N = 924;
CRLB_a1 = (1 - a1^2) ./ Ngrid;

figure;
imagesc(N_vals, sigma2_vals, CRLB_sigma2);
set(gca, 'YDir', 'normal', 'XScale', 'log', 'YScale', 'log', 'ColorScale', 'log');
colormap(turbo(256));
colorbar;
xlabel('N', 'FontSize', 14);
ylabel('\sigma^2', 'FontSize', 14);
title('CRLB for \sigma^2', 'FontSize', 14);
grid on;

figure;
imagesc(N_vals, sigma2_vals, CRLB_a1);
set(gca, 'YDir', 'normal', 'XScale', 'log', 'YScale', 'log', 'ColorScale', 'log');
colormap(turbo(256));
colorbar;
xlabel('N', 'FontSize', 14);
ylabel('\sigma^2', 'FontSize', 14);
title(sprintf('CRLB for a1 = %.5f', a1), 'FontSize', 14);
grid on;
