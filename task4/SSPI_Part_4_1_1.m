clear; clc; close all;

N = 1000;
sigma = 0.1;
b = [1, 2, 3, 2, 1];
a = 1;
Nw = numel(b) - 1; % Wiener filter order

rng(42); 

x = randn(N, 1);
y_raw = filter(b, a, x);

% Normalise y so that var(y) = 1 (population variance)
scale = 1 / std(y_raw, 1);
y = scale * y_raw;

eta = sigma * randn(N, 1);
z = y + eta;

% Compute SNR
signalPower = mean(y .^ 2);
noisePower = mean(eta .^ 2);
snrDb = 10 * log10(signalPower / noisePower);

fprintf('N = %d\n', N);
fprintf('var(y) after normalisation: %.6f\n', var(y, 1));
fprintf('var(eta): %.6f (target %.6f)\n', var(eta, 1), sigma^2);
fprintf('SNR(z) = %.3f dB\n\n', snrDb);
fprintf('Output scale c = %.6f, so ideal target coeffs are c*b (not b).\n\n', scale);

lagsRxx = -Nw:Nw;
rxx = corr_est(x, x, lagsRxx);

zeroLagIdx = find(lagsRxx == 0);
rxx_pos = rxx(zeroLagIdx:end);
rxx_neg = rxx(zeroLagIdx:-1:1);  
Rxx = toeplitz(rxx_pos, rxx_neg);

lagsRzx = -Nw:0;
rzx = corr_est(z, x, lagsRzx);

pzx = rzx(end:-1:1).';

disp('Rxx =');
disp(Rxx);
disp('pzx =');
disp(pzx);
fprintf('\n');

% Wiener solution
wopt = Rxx \ pzx;

% Unknown-system coefficients 
b_norm = (scale * b).';
relErr = norm(wopt - b_norm) / norm(b_norm);

disp('Optimal Wiener coefficients wopt^T =');
disp(wopt.');
disp('Normalised unknown-system coefficients b_norm^T =');
disp(b_norm.');
fprintf('Relative L2 error = %.6f\n', relErr);

if relErr < 0.20
    fprintf('Conclusion: wopt is similar to the unknown system coefficients.\n');
else
    fprintf('Conclusion: wopt is not very close (finite-sample estimation error).\n');
end

figure;
stem(0:Nw, b_norm, 'filled', 'DisplayName', 'Unknown system (normalised)');
hold on;
stem(0:Nw, wopt, 'r', 'DisplayName', 'Wiener estimate');
hold off;
grid on;
xlabel('Coefficient index');
ylabel('Coefficient value');
title('Unknown system vs Wiener-optimal estimate');
legend('Location', 'best');

% Noise power sweep
sigma2_vals = logspace(-1, 1, 6);
snrDb_vals = zeros(size(sigma2_vals));
relErr_vals = zeros(size(sigma2_vals));
fprintf('\nNoise power sweep (Nw=%d):\n', Nw);
fprintf('sigma^2\t\tSNR(dB)\t\trelErr\n');
for i = 1:numel(sigma2_vals)
    sigma_i = sqrt(sigma2_vals(i));
    eta_i = sigma_i * randn(N, 1);
    z_i = y + eta_i;

    signalPower_i = mean(y .^ 2);
    noisePower_i = mean(eta_i .^ 2);
    snrDb_i = 10 * log10(signalPower_i / noisePower_i);

    rzx_i = corr_est(z_i, x, lagsRzx);
    pzx_i = rzx_i(end:-1:1).';
    wopt_i = Rxx \ pzx_i;
    relErr_i = norm(wopt_i - b_norm) / norm(b_norm);

    fprintf('%.4f\t\t%.3f\t\t%.6f\n', sigma2_vals(i), snrDb_i, relErr_i);
    snrDb_vals(i) = snrDb_i;
    relErr_vals(i) = relErr_i;
end

figure;
semilogx(sigma2_vals, snrDb_vals, 'o-', 'LineWidth', 1.2);
grid on;
xlabel('\sigma^2');
ylabel('SNR (dB)');
title('SNR vs Noise Variance');

figure;
semilogx(sigma2_vals, relErr_vals, 's-', 'LineWidth', 1.2);
grid on;
xlabel('\sigma^2');
ylabel('Relative L2 error');
title('Wiener Error vs Noise Variance');

Nw_list = [4, 8, 16, 32];
fprintf('\nEffect of larger Nw (sigma=%.2f):\n', sigma);
fprintf('Nw\t\trelErr\n');
for i = 1:numel(Nw_list)
    Nw_i = Nw_list(i);
    lagsRxx_i = -Nw_i:Nw_i;
    rxx_i = corr_est(x, x, lagsRxx_i);
    zeroLagIdx_i = find(lagsRxx_i == 0);
    rxx_pos_i = rxx_i(zeroLagIdx_i:end);
    rxx_neg_i = rxx_i(zeroLagIdx_i:-1:1);
    Rxx_i = toeplitz(rxx_pos_i, rxx_neg_i);

    lagsRzx_i = -Nw_i:0;
    rzx_i = corr_est(z, x, lagsRzx_i);
    pzx_i = rzx_i(end:-1:1).';

    wopt_i = Rxx_i \ pzx_i;
    b_pad = [b_norm; zeros(Nw_i + 1 - numel(b_norm), 1)];
    relErr_i = norm(wopt_i - b_pad) / max(norm(b_pad), eps);

    fprintf('%d\t\t%.6f\n', Nw_i, relErr_i);
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
