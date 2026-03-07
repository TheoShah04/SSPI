clear; clc;
S = load('sunspot.dat');
if isstruct(S)
    fn = fieldnames(S);
    D = S.(fn{1});
else
    D = S;
end

if size(D, 2) >= 2
    x = D(:, 2);
else
    x = D(:);
end

x = x(:);
x_mean = mean(x);
x = x - x_mean; % mean-center
maxP = 10;
Mfactor = 4; % overdetermined LS: M = 4p (M >= p)

p_out = (1:maxP).';
M_out = zeros(maxP, 1);
diffNorm = zeros(maxP, 1); % ||a_LS - a_YW||_2 (reference)
a_ls_all = cell(maxP, 1);
a_yw_all = cell(maxP, 1);

for p = 1:maxP
    M = Mfactor * p;
    M_out(p) = M;

    % Need ACF values for lags 0..M
    r = xcorr(x, M, 'biased');
    rpos = r(M+1:end);          % r(0), r(1), ..., r(M)

    s = rpos(2:M+1);            % [r(1) ... r(M)]^T
    H = zeros(M, p);
    for k = 1:M
        for i = 1:p
            lag = k - i;
            H(k, i) = rpos(abs(lag) + 1); % r(-m) = r(m) for real processes
        end
    end

    a_ls = H \ s;
    a_ls_all{p} = a_ls;

    % Yule-Walker reference from aryule
    A = aryule(x, p);  
    a_yw = -A(2:end).';
    a_yw_all{p} = a_yw;

    diffNorm(p) = norm(a_ls - a_yw, 2);

    fprintf('p = %2d, M = %2d\n', p, M);
    fprintf('  a_LS = [%s]\n', sprintf(' %.5f', a_ls));
    fprintf('  a_YW = [%s]\n', sprintf(' %.5f', a_yw));
    fprintf('  ||a_LS - a_YW||_2 = %.5e\n\n', diffNorm(p));
end

T = table(p_out, M_out, diffNorm, ...
    'VariableNames', {'p', 'M', 'norm_aLS_minus_aYW'});
disp(T);

figure;
plot(p_out, diffNorm, '-o', 'LineWidth', 1.2);
grid on;
xlabel('AR order p', 'FontSize', 14);
ylabel('||a_{LS} - a_{YW}||_2', 'FontSize', 14);
title('Difference Between LS and Yule-Walker Estimates', 'FontSize', 14);
set(gca, 'FontSize', 14);

coeffMat = NaN(maxP, maxP); % row = p, col = coefficient index i
for p = 1:maxP
    coeffMat(p, 1:p) = a_ls_all{p}.';
end

maxAbsCoeff = max(abs(coeffMat(:)), [], 'omitnan');

figure;
hold on;
clr = turbo(maxP);
plotOrders = [1 2 4 8 10];
for p = plotOrders
    hp = stem(1:p, a_ls_all{p}, 'LineWidth', 1.0, 'DisplayName', sprintf('p=%d', p));
    hp.Color = clr(p, :);
    hp.MarkerFaceColor = clr(p, :);
    hp.MarkerSize = 4;
end
hold off;
grid on;
xlim([0.5 maxP + 0.5]);
ylim([-1.05 * maxAbsCoeff, 1.05 * maxAbsCoeff]);
xlabel('Coefficient index i', 'FontSize', 14);
ylabel('a_i', 'FontSize', 14);
title('LS AR Coefficients (Overlay Stem Plot)', 'FontSize', 14);
legend('show', 'Location', 'eastoutside');
set(gca, 'FontSize', 14);

N = numel(x);
rmse = NaN(maxP, 1);
mse = NaN(maxP, 1);
errCell = cell(maxP, 1);

for p = 1:maxP
    a = a_ls_all{p}(:);
    y = x(p+1:end);
    Xreg = zeros(N - p, p);
    for i = 1:p
        Xreg(:, i) = x(p+1-i:N-i);
    end

    yhat = Xreg * a;
    e = y - yhat;
    errCell{p} = e;
    mse(p) = mean(e.^2);
    rmse(p) = sqrt(mse(p));
end

[~, bestP] = min(rmse);
fprintf('\nBest model order by minimum RMSE: p = %d (RMSE = %.5f)\n', bestP, rmse(bestP));

figure;
plot(1:maxP, rmse, '-o', 'LineWidth', 1.3);
grid on;
xlabel('Model order p', 'FontSize', 14);
ylabel('Approximation RMSE', 'FontSize', 14);
title('AR Model Approximation Error vs Order', 'FontSize', 14);
set(gca, 'FontSize', 14);

figure;
hold on;
plotOrdersErr = [1 2 4 8 10];
for p = plotOrdersErr
    plot(errCell{p}, 'LineWidth', 1.0, 'DisplayName', sprintf('p=%d', p));
end
hold off;
grid on;
xlabel('Sample index n', 'FontSize', 14);
ylabel('Approximation error e[n]', 'FontSize', 14);
title('Approximation Error to Original Data for Selected AR Orders', 'FontSize', 14);
legend('show', 'Location', 'eastoutside');
set(gca, 'FontSize', 14);

% Power spectra
Nfft = 1024;
clrPSD = turbo(maxP); % force unique, order-consistent colors

figure;
hold on;
for p = 1:maxP
    a = a_ls_all{p}(:);
    den = [1, -a.']; % A(z) = 1 - sum_i a_i z^{-i}
    num = sqrt(mse(p));
    [h, w] = freqz(num, den, Nfft);
    semilogy(w/(2*pi), abs(h).^2 + eps, 'LineWidth', 1.0, ...
        'Color', clrPSD(p, :), 'DisplayName', sprintf('AR(%d)', p));
end
hold off;
grid on;
xlim([0 0.5]);
xlabel('Normalised frequency', 'FontSize', 16);
ylabel('PSD', 'FontSize', 16);
title('Sunspot AR(p) Model-Based Power Spectra (LS)', 'FontSize', 16);
legend('show', 'Location', 'northeast');
set(gca, 'FontSize', 14);

figure;
hold on;
plotOrdersPSD = [1 2 4 8 10];
for p = plotOrdersPSD
    a = a_ls_all{p}(:);
    den = [1, -a.'];
    num = sqrt(mse(p));
    [h, w] = freqz(num, den, Nfft);
    semilogy(w/(2*pi), abs(h).^2 + eps, 'LineWidth', 1.4, ...
        'Color', clrPSD(p, :), 'DisplayName', sprintf('AR(%d)', p));
end
hold off;
grid on;
xlim([0 0.5]);
xlabel('Normalized frequency', 'FontSize', 14);
ylabel('PSD', 'FontSize', 14);
title('Sunspot AR(p) Power Spectra for Selected Orders', 'FontSize', 14);
legend('show', 'Location', 'eastoutside');
set(gca, 'FontSize', 14);
