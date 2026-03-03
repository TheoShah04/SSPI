clear; clc;

load sunspot.dat
x = sunspot(:, 2);
x = x(:);

xz = (x - mean(x)) / std(x);

pVals = [1 2 10];
mVals = [1 2 5 10];
N = length(xz);

% Store error statistics for each (p, m)
MSE = zeros(numel(pVals), numel(mVals));
Bias2 = zeros(numel(pVals), numel(mVals));
VarE = zeros(numel(pVals), numel(mVals));

% Precompute AR coefficients
aVals = cell(numel(pVals), 1);
for ip = 1:numel(pVals)
    aVals{ip} = ar_yw(xz, pVals(ip));
end

figure;
tiledlayout(2, 2, 'TileSpacing', 'compact');

maxp = max(pVals);
colors = lines(numel(pVals));

for im = 1:numel(mVals)
    m = mVals(im);
    idx = (maxp + m):N;

    nexttile;
    plot(idx, xz(idx), 'k-', 'LineWidth', 1.0, 'DisplayName', 'Actual');
    hold on;

    for ip = 1:numel(pVals)
        p = pVals(ip);
        pred = mstep_predict(xz, aVals{ip}, p, m);
        e = xz(idx) - pred(idx);
        MSE(ip, im) = mean(e.^2);
        Bias2(ip, im) = mean(e)^2;
        VarE(ip, im) = var(e, 1);
        plot(idx, pred(idx), '-', 'LineWidth', 1.5, 'Color', colors(ip, :), ...
            'DisplayName', sprintf('AR(%d)', p));
    end

    hold off;
    xlabel('Sample number', 'FontSize', 20);
    ylabel('Standardized value', 'FontSize', 16);
    title(sprintf('m = %d', m), 'FontSize', 20);
    grid on;
    if im == 1
        legend('show', 'Location', 'best');
    end
end

disp('MSE for each model order p and horizon m (rows: p, columns: m)');
disp(array2table(MSE, 'VariableNames', compose('m_%d', mVals), ...
    'RowNames', compose('p_%d', pVals)));

% Bias-variance tradeoff plots
figure;
tiledlayout(1, 2, 'TileSpacing', 'compact');

nexttile;
hold on;
for ip = 1:numel(pVals)
    plot(mVals, Bias2(ip, :), '-o', 'LineWidth', 1.5, ...
        'Color', colors(ip, :), 'DisplayName', sprintf('AR(%d)', pVals(ip)));
end
hold off;
xlabel('Prediction horizon m', 'FontSize', 16);
ylabel('Bias^2', 'FontSize', 16);
title('Bias^2 vs. horizon', 'FontSize', 16);
grid on;
legend('show', 'Location', 'best');

nexttile;
hold on;
for ip = 1:numel(pVals)
    plot(mVals, VarE(ip, :), '-o', 'LineWidth', 1.5, ...
        'Color', colors(ip, :), 'DisplayName', sprintf('AR(%d)', pVals(ip)));
end
hold off;
xlabel('Prediction horizon m', 'FontSize', 16);
ylabel('Variance of error', 'FontSize', 16);
title('Variance vs. horizon', 'FontSize', 16);
grid on;
legend('show', 'Location', 'best');

function a = ar_yw(x, p)
    % Yule-Walker AR(p) coefficients for x[n] = a1 x[n-1] + ... + ap x[n-p] + e[n]
    rxx = xcorr(x, p, 'unbiased');
    r0p = rxx(p+1:end);           % lags 0..p
    R = toeplitz(r0p(1:p));       % r0..r_{p-1}
    rvec = r0p(2:p+1);            % r1..rp
    a = R \ rvec;
end

function pred = mstep_predict(x, a, p, m)
    % m-step ahead prediction using recursive AR updates
    N = length(x);
    pred = nan(N, 1);
    for n = p:N-m
        state = x(n:-1:n-p+1);
        xpred = 0;
        for k = 1:m
            xpred = a' * state;
            state = [xpred; state(1:end-1)];
        end
        pred(n + m) = xpred;
    end
end
