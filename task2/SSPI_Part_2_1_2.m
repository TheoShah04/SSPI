N = 1000;
rng(42);
x = randn(N, 1);

% Base MA(9) with unit coefficients
b = ones(9, 1);
a = 1;

y = filter(b, a, x);
[r, lags] = xcorr(y, 'unbiased');

idx = (lags >= -20) & (lags <= 20);
figure;
stem(lags(idx), r(idx), 'filled');
xlabel('Lag \tau');
ylabel('Unbiased ACF estimate');
title('ACF of MA(9) filtered WGN');
grid on;

orders = [1 3 5 9 15];
figure;
hold on;
for k = 1:numel(orders)
    M = orders(k);
    b = ones(M, 1);
    y = filter(b, 1, x);
    [r, lags] = xcorr(y, 'unbiased');
    idx = (lags >= -20) & (lags <= 20);
    stem(lags(idx), r(idx), 'filled', 'DisplayName', sprintf('order=%d', M-1));
end
hold off;
xlabel('Lag \tau', 'FontSize', 14);
ylabel('Unbiased ACF estimate', 'FontSize', 14);
title('ACF estimate vs MA Order (unit coefficients)');
legend('show', 'Location', 'northeast');
grid on;

M = 9;
coeffs = {
    0.5 * ones(M,1), ...
    1.0 * ones(M,1), ...      
    2.0 * ones(M,1) ...       
};
labels = {
    'Rectangular (0.5)',
    'Rectangular (1.0)',
    'Rectangular (2.0)'
};

figure;
hold on;
for k = 1:numel(coeffs)
    b = coeffs{k};
    y = filter(b, 1, x);
    [r, lags] = xcorr(y, 'unbiased');
    idx = (lags >= -20) & (lags <= 20);
    stem(lags(idx), r(idx), 'filled', 'DisplayName', labels{k});
end
hold off;
xlabel('Lag \tau', 'FontSize', 14);
ylabel('Unbiased ACF estimate', 'FontSize', 14);
title('ACF estimate vs Coefficient Value (Order 9)');
legend('show', 'Location', 'northeast');
grid on;
