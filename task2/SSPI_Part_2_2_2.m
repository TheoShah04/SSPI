N = 5000;
rng(42);
x = randn(N, 1);           % WGN input
sigma2 = var(x, 1);        % population variance estimate

orders = [2 4 8 12];        % filter orders (order = M-1)

% Overlap all orders on the same figure
figure;
hold on;
for k = 1:numel(orders)
    p = orders(k);          % order
    M = p + 1;              % length
    h = ones(M, 1);         % MA coefficients

    y = filter(h, 1, x);

    % Unbiased CCF estimate
    [r_xy, lags] = xcorr(x, y, 'unbiased');

    % Theoretical CCF: sigma^2 * h(-tau)
    r_th = zeros(size(r_xy));
    neg = (lags <= 0) & (lags >= -(M-1));
    r_th(neg) = sigma2 * h(-lags(neg) + 1);

    % Plot around zero lag
    idx = (lags >= -20) & (lags <= 20);
    stem(lags(idx), r_xy(idx), 'filled', ...
        'DisplayName', sprintf('Est (order=%d)', p));

    % Reconstruct h from CCF (system identification idea)
    h_hat = flipud(r_xy(neg) / sigma2);   % map negative lags to n=0..p
    rel_err = norm(h_hat - h) / norm(h);
    fprintf('Order %d: relative h-hat error = %.3f\n', p, rel_err);
end
hold off;
xlabel('Lag \tau', 'FontSize', 14);
ylabel('CCF', 'FontSize', 14);
title('CCF for varying filter orders', 'FontSize', 14);
legend('show', 'Location', 'northeast');
grid on;

% Summary:
% Increasing order p increases the CCF support length (from -p to 0),
% i.e., a wider rectangular block for unit coefficients.
