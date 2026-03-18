clear; clc;

S = load('NASDAQ.mat');
if isfield(S, 'NASDAQ')
    NASDAQ = S.NASDAQ;
else
    fn = fieldnames(S);
    NASDAQ = S.(fn{1});
end

close = NASDAQ.Close; 
dates = NASDAQ.Date;

close = close(:);
dates = dates(:);

disp('First 10 rows of NASDAQ data (Date, Close):');
disp(table(dates(1:10), close(1:10), 'VariableNames', {'Date','Close'}));
fprintf('Number of samples (raw close prices): %d\n', numel(close));

% Time series plot of NASDAQ close
figure;
plot(dates, close, 'LineWidth', 1);
xlabel('Date', 'FontSize', 14);
ylabel('Close', 'FontSize', 14);
title('NASDAQ Closing Prices', 'FontSize', 14);
grid on;

% Daily return representations
r_simple = 100 * (close(2:end) - close(1:end-1)) ./ close(1:end-1);
r_log = 100 * log(close(2:end) ./ close(1:end-1));
r_abs = abs(r_simple);

returns = {r_simple, r_log, r_abs};
return_labels = {'Simple % return', 'Log % return', 'Absolute % return'};

for i = 1:numel(returns)
    ri = returns{i};
    ri = ri(~isnan(ri));
    ri = ri - mean(ri);
    returns{i} = ri;
end

r = returns{1};
N = length(r);

% ACF plot for simple returns
maxp = 10;
figure;
autocorr(r, 'NumLags', maxp);
xlabel('Lag', 'FontSize', 14);
ylabel('ACF', 'FontSize', 14);
title('Autocorrelation of NASDAQ Returns (Simple %)', 'FontSize', 14);
grid on;

% PACF plot for closing prices
figure;
parcorr(close - mean(close), 'NumLags', maxp);
xlabel('Lag', 'FontSize', 14);
ylabel('PACF', 'FontSize', 14);
title('PACF of NASDAQ Closing Prices', 'FontSize', 14);
grid on;

% PACF plots for each representation
for i = 1:numel(returns)
    figure;
    parcorr(returns{i}, 'NumLags', maxp);
    xlabel('Lag', 'FontSize', 14);
    ylabel('PACF', 'FontSize', 14);
    title(sprintf('PACF of NASDAQ Returns - %s', return_labels{i}), 'FontSize', 14);
    grid on;
end

p_list = (1:maxp)';
AIC = NaN(size(p_list));
MDL = NaN(size(p_list));
a1 = NaN;

for i = 1:numel(p_list)
    p = p_list(i);
    mdl = arima(p, 0, 0);
    mdl.Constant = 0;
    mdl.Variance = NaN;
    try
        [EstMdl, ~, logL] = estimate(mdl, r, 'Display', 'off');
        numParams = p + 1; % AR coefficients + variance
        [AIC(i), BIC] = aicbic(logL, numParams, N);
        MDL(i) = BIC;
        if p == 1
            a1 = EstMdl.AR{1};
        end
    catch ME
        warning('AR(%d) estimation failed: %s', p, ME.message);
    end
end

[~, idx_aic] = min(AIC);
[~, idx_mdl] = min(MDL);
p_aic = p_list(idx_aic);
p_mdl = p_list(idx_mdl);

figure;
plot(1:maxp, AIC, '-o', 'DisplayName', 'AIC');
hold on;
plot(1:maxp, MDL, '-s', 'DisplayName', 'MDL');
hold off;
xlabel('Model order p', 'FontSize', 14);
ylabel('Criterion value', 'FontSize', 14);
title('AIC and MDL for AR(p) Models', 'FontSize', 14);
legend('show', 'Location', 'northeast');
grid on;

fprintf('Estimated a1 (AR(1)) = %.6f\n', a1);
