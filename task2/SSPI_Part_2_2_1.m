N = 1000;
rng(42);
x = randn(N, 1);           
b = ones(11, 1);       

y = filter(b, 1, x);  

% Unbiased CCF estimate
[r_xy, lags] = xcorr(x, y, 'unbiased');

idx = (lags >= -20) & (lags <= 20);
figure;
stem(lags(idx), r_xy(idx), 'filled');
xlabel('Lag \tau', 'FontSize', 14);
ylabel('Unbiased CCF estimate', 'FontSize', 14);
title('CCF between x and y', 'FontSize', 14);
grid on;

