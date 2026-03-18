N = 5000;
rng(42);
x = randn(N, 1);         
sigma2 = var(x, 1);      

orders = [2 4 8 12];      

figure;
hold on;
for k = 1:numel(orders)
    p = orders(k);      
    M = p + 1;           
    h = ones(M, 1); % MA coefficients

    y = filter(h, 1, x);

    % Unbiased CCF estimate
    [r_xy, lags] = xcorr(x, y, 'unbiased');

    % Theoretical CCF
    r_th = zeros(size(r_xy));
    neg = (lags <= 0) & (lags >= -(M-1));
    r_th(neg) = sigma2 * h(-lags(neg) + 1);

    idx = (lags >= -20) & (lags <= 20);
    stem(lags(idx), r_xy(idx), 'filled', ...
        'DisplayName', sprintf('Est (order=%d)', p));

    h_hat = flipud(r_xy(neg) / sigma2);
    rel_err = norm(h_hat - h) / norm(h);
    fprintf('Order %d: relative h-hat error = %.3f\n', p, rel_err);
end
hold off;
xlabel('Lag \tau', 'FontSize', 14);
ylabel('CCF', 'FontSize', 14);
title('CCF for varying filter orders', 'FontSize', 14);
legend('show', 'Location', 'northeast');
grid on;