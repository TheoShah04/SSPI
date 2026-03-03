N = 100000;
rng(42);
x = randn(N, 1);

[r, lags] = xcorr(x, 'unbiased');

figure;
stem(lags, r, 'filled');
axis([-99999 99999 (min(r)-0.5) (max(r)+0.5)]);
xlabel('Lag \tau');
ylabel('Unbiased ACF estimate');
title('Unbiased ACF estimate for WGN (N=10,000)');
grid on;