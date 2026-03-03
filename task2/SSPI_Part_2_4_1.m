% AR(1) sufficiency for NASDAQ daily returns using PACF and AIC/MDL
clear; clc;

% Load NASDAQ data
S = load('NASDAQ.mat');
if isfield(S, 'NASDAQ')
    NASDAQ = S.NASDAQ;
else
    fn = fieldnames(S);
    NASDAQ = S.(fn{1});
end

% Extract closing prices and dates
if istable(NASDAQ)
    if any(strcmp('Close', NASDAQ.Properties.VariableNames))
        close = NASDAQ.Close;
    else
        close = NASDAQ{:, end};
    end
    if any(strcmp('Date', NASDAQ.Properties.VariableNames))
        dates = NASDAQ.Date;
    else
        dates = (1:height(NASDAQ))';
    end
elseif isstruct(NASDAQ)
    if isfield(NASDAQ, 'Close')
        close = NASDAQ.Close;
    else
        fn = fieldnames(NASDAQ);
        close = NASDAQ.(fn{end});
    end
    if isfield(NASDAQ, 'Date')
        dates = NASDAQ.Date;
    else
        dates = (1:numel(close))';
    end
else
    close = NASDAQ(:, end);
    dates = NASDAQ(:, 1);
end

close = close(:);
dates = dates(:);

% Snippet of raw data
disp('First 10 rows of NASDAQ data (Date, Close):');
disp(table(dates(1:10), close(1:10), 'VariableNames', {'Date','Close'}));

% Daily log-returns: log(p_n / p_{n-1})
r = log(close(2:end) ./ close(1:end-1));
r = r(~isnan(r));
N = length(r);

% Zero-mean returns
r = r - mean(r);

% PACF via Yule-Walker for orders 1..maxp
maxp = 10;
pacf = zeros(maxp, 1);
sigma2 = zeros(maxp + 1, 1); % include p = 0
a1 = NaN;

% p = 0 (white noise) variance
sigma2(1) = mean(r.^2);

for p = 1:maxp
    rxx = xcorr(r, p, 'unbiased');          % lags -p..p
    r0p = rxx(p+1:end);                     % lags 0..p
    R = toeplitz(r0p(1:p));                 % r0..r_{p-1}
    rvec = r0p(2:p+1);                      % r1..rp
    a = R \ rvec;                           % AR(p) coefficients
    pacf(p) = a(end);                       % last coefficient = PACF(p)
    sigma2(p + 1) = r0p(1) - a' * rvec;     % prediction error variance
    if p == 1
        a1 = a(1);
    end
end

% Information criteria (standard forms)
p_list = (0:maxp)';
AIC = log(sigma2) + (2 * p_list) / N;
MDL = log(sigma2) + (p_list * log(N)) / N;
[~, idx_aic] = min(AIC);
[~, idx_mdl] = min(MDL);
p_aic = p_list(idx_aic);
p_mdl = p_list(idx_mdl);

% PACF plot with 95% bounds
figure;
stem(1:maxp, pacf, 'filled');
hold on;
conf = 1.96 / sqrt(N);
yline(conf, 'r--');
yline(-conf, 'r--');
hold off;
xlabel('Lag', 'FontSize', 14);
ylabel('PACF', 'FontSize', 14);
title('Partial Autocorrelation of NASDAQ Returns', 'FontSize', 14);
grid on;

figure;
plot(0:maxp, AIC, '-o', 'DisplayName', 'AIC');
hold on;
plot(0:maxp, MDL, '-s', 'DisplayName', 'MDL');
hold off;
xlabel('Model order p', 'FontSize', 14);
ylabel('Criterion value', 'FontSize', 14);
title('AIC and MDL for AR(p) Models', 'FontSize', 14);
legend('show', 'Location', 'northeast');
grid on;

fprintf('Estimated a1 (AR(1)) = %.6f\n', a1);
