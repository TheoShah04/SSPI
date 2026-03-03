clear; clc;

N_vals = 1:50:1001;
sigma2_vals = 1:50:1001;

[Ngrid, Sgrid] = meshgrid(N_vals, sigma2_vals);

% CRLB for sigma^2: var(sigma^2_hat) >= 2*sigma^4 / N
% Here Sgrid represents sigma^2, so sigma^4 = (sigma^2)^2
CRLB_sigma2 = 2 .* (Sgrid .^ 2) ./ Ngrid;

% CRLB for a1: var(a1_hat) >= (1 - a1^2) / N
% Set a1 to your estimated AR(1) coefficient if available
a1 = -0.005581;
CRLB_a1 = (1 - a1^2) ./ Ngrid;

figure;
h1 = heatmap(N_vals, sigma2_vals, CRLB_sigma2);
h1.XLabel = 'N';
h1.YLabel = '\sigma^2';
h1.Title = 'CRLB for \sigma^2';
h1.FontSize = 14;

figure;
h2 = heatmap(N_vals, sigma2_vals, CRLB_a1);
h2.XLabel = 'N';
h2.YLabel = '\sigma^2';
h2.Title = sprintf('CRLB for a1 = %.2f', a1);
h2.FontSize = 14;
