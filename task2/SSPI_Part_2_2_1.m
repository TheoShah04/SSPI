% CCF estimate for x and y from Section 2.1
N = 1000;
rng(42);
x = randn(N, 1);            % WGN
b = ones(11, 1);             % MA(9) coefficients

y = filter(b, 1, x);        % Filtered output

% Unbiased CCF estimate
[r_xy, lags] = xcorr(x, y, 'unbiased');

% Plot lags in [-20, 20]
idx = (lags >= -20) & (lags <= 20);
figure;
stem(lags(idx), r_xy(idx), 'filled');
xlabel('Lag \tau', 'FontSize', 14);
ylabel('Unbiased CCF estimate', 'FontSize', 14);
title('CCF between x and y', 'FontSize', 14);
grid on;

