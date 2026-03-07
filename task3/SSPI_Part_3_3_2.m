clear; clc;
S = load('sunspot.dat');
if isstruct(S)
    fn = fieldnames(S);
    D = S.(fn{1});
else
    D = S;
end

if size(D, 2) >= 2
    x_all = D(:, 2);
else
    x_all = D(:);
end
x_all = x_all(:);

N_vals = (10:5:250).';
maxP = 10;
Mfactor = 5;                  % use M = 4p when feasible
N_ref = 100;
N_vals(end);          % choose optimal p at largest N in range

if numel(x_all) < N_ref
    error('sunspot.dat has fewer than %d samples.', N_ref);
end

% --- Step 1: find optimal order p* at N_ref by LS prediction MSE ---
x_ref = x_all(1:N_ref);
x_ref = x_ref - mean(x_ref);

mse_by_p = NaN(maxP, 1);
J_by_p = NaN(maxP, 1);

for p = 1:maxP
    if N_ref <= p
        continue;
    end

    M = min(Mfactor * p, N_ref - 1);
    if M < p
        continue;
    end

    r = xcorr(x_ref, M, 'biased');
    rpos = r(M + 1:end);              % lags 0..M

    svec = rpos(2:M+1);               % [r(1) ... r(M)]^T
    H = zeros(M, p);
    for k = 1:M
        for i = 1:p
            lag = k - i;
            H(k, i) = rpos(abs(lag) + 1); % r(-m)=r(m)
        end
    end

    a_ls = H \ svec;                  % LS solution minimizing J
    J_by_p(p) = norm(svec - H * a_ls)^2;

    y = x_ref(p+1:end);
    Xreg = zeros(N_ref - p, p);
    for i = 1:p
        Xreg(:, i) = x_ref(p+1-i:N_ref-i);
    end
    e = y - Xreg * a_ls;
    mse_by_p(p) = mean(e.^2);
end

[~, p_opt] = min(mse_by_p);
fprintf('Optimal AR order at N=%d: p*=%d (MSE=%.6f)\n', N_ref, p_opt, mse_by_p(p_opt));

% Order-selection plot (shows why p_opt is chosen)
figure;
plot(1:maxP, mse_by_p, '-o', 'LineWidth', 1.3, 'DisplayName', 'MSE vs order');
hold on;
plot(p_opt, mse_by_p(p_opt), 'rp', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'r', 'DisplayName', sprintf('Selected p=%d', p_opt));
hold off;
grid on;
xlabel('AR model order p', 'FontSize', 14);
ylabel('Validation MSE at N=250', 'FontSize', 14);
title('AR Order Selection by Minimum MSE', 'FontSize', 14);
legend('show', 'Location', 'northeast');
set(gca, 'FontSize', 14);

% --- Step 2: approximation error vs N for fixed p* ---
mse_vs_N = NaN(size(N_vals));
J_vs_N = NaN(size(N_vals));
M_vs_N = NaN(size(N_vals));

for idx = 1:numel(N_vals)
    N = N_vals(idx);
    if N <= p_opt
        continue;
    end

    xN = x_all(1:N);
    xN = xN - mean(xN);

    M = min(Mfactor * p_opt, N - 1);
    if M < p_opt
        continue;
    end
    M_vs_N(idx) = M;

    r = xcorr(xN, M, 'biased');
    rpos = r(M + 1:end);              % lags 0..M

    svec = rpos(2:M+1);
    H = zeros(M, p_opt);
    for k = 1:M
        for i = 1:p_opt
            lag = k - i;
            H(k, i) = rpos(abs(lag) + 1);
        end
    end

    a_ls = H \ svec;
    J_vs_N(idx) = norm(svec - H * a_ls)^2;

    y = xN(p_opt+1:end);
    Xreg = zeros(N - p_opt, p_opt);
    for i = 1:p_opt
        Xreg(:, i) = xN(p_opt+1-i:N-i);
    end
    e = y - Xreg * a_ls;
    mse_vs_N(idx) = mean(e.^2);
end

% Suggest practical N: first N within 5% of minimum MSE
valid = isfinite(mse_vs_N);
mse_min = min(mse_vs_N(valid));
idx_optN = find(valid & (mse_vs_N <= 1.05 * mse_min), 1, 'first');
N_opt = N_vals(idx_optN);
fprintf('Suggested data length (within 5%% of min MSE): N=%d\n', N_opt);
fprintf('Minimum MSE in range: %.6f\n', mse_min);

% --- Plots ---
figure;
plot(N_vals, mse_vs_N, '-o', 'LineWidth', 1.3);
grid on;
xlabel('Data length N', 'FontSize', 14);
ylabel('Approximation error (MSE)', 'FontSize', 14);
title(sprintf('MSE vs N for Optimal AR(%d) Model', p_opt), 'FontSize', 14);
set(gca, 'FontSize', 14);

figure;
plot(N_vals, J_vs_N, '-s', 'LineWidth', 1.3);
grid on;
xlabel('Data length N', 'FontSize', 14);
ylabel('Minimum LS cost J', 'FontSize', 14);
title(sprintf('LS Cost J vs N for AR(%d)', p_opt), 'FontSize', 14);
set(gca, 'FontSize', 14);
