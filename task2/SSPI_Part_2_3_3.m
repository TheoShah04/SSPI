% Sunspot PACF via Yule-Walker up to p = 10, with and without standardization
clear; clc;

load sunspot.dat
x = sunspot(:, 2);
x = x(:);

maxp = 10;

x0 = x; 

% Standardised
xz = (x - mean(x)) / std(x);

[pacf0, conf0] = pacf_yw(x0, maxp);
[pacfz, confz] = pacf_yw(xz, maxp);

% (95% bound)
p0 = last_significant_lag(pacf0, conf0);
pz = last_significant_lag(pacfz, confz);

% Plots
figure;
subplot(2,1,1);
stem(1:maxp, pacf0, 'filled');
hold on;
yline(conf0, 'r--');
yline(-conf0, 'r--');
hold off;
xlabel('Lag');
ylabel('PACF');
title('Sunspot PACF (non-standardised)');
grid on;

subplot(2,1,2);
stem(1:maxp, pacfz, 'filled');
hold on;
yline(confz, 'r--');
yline(-confz, 'r--');
hold off;
xlabel('Lag');
ylabel('PACF');
title('Sunspot PACF (zero-mean, unit variance)');
grid on;

fprintf('Suggested order (zero-mean): p = %d\n', p0);
fprintf('Suggested order (standardized): p = %d\n', pz);

if p0 == pz
    fprintf('Standardization does not change the PACF shape; model order is unchanged.\n');
else
    fprintf('Differences are due to mean/scale effects or finite-sample bias.\n');
end

function [pacf, conf] = pacf_yw(x, maxp)
    x = x(:);
    N = length(x);
    pacf = zeros(maxp, 1);

    for p = 1:maxp
        rxx = xcorr(x, p, 'unbiased');      % lags -p..p
        r0p = rxx(p+1:end);                 % lags 0..p
        R = toeplitz(r0p(1:p));             % r0..r_{p-1}
        rvec = r0p(2:p+1);                  % r1..rp
        a = R \ rvec;                       % Yule-Walker AR(p)
        pacf(p) = a(end);                   % last coefficient = PACF(p)
    end

    conf = 1.96 / sqrt(N);
end

function p = last_significant_lag(pacf, conf)
    idx = find(abs(pacf) > conf);
    if isempty(idx)
        p = 0;
    else
        p = idx(end);
    end
end
