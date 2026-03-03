clear; clc;

load sunspot.dat
N = size(sunspot,1);
fprintf('Size of dataset %d\n', N);
x = sunspot(:, 2);
x = x(:);

xz = (x - mean(x)) / std(x);

maxp = 10;
N = length(xz);

MDL = zeros(maxp, 1);
AIC = zeros(maxp, 1);
AICc = zeros(maxp, 1);
Ep = zeros(maxp, 1);

for p = 1:maxp
    % Yule-Walker AR(p)
    rxx = xcorr(xz, p, 'unbiased');      % lags -p..p
    r0p = rxx(p+1:end);                  % lags 0..p
    R = toeplitz(r0p(1:p));              % r0..r_{p-1}
    rvec = r0p(2:p+1);                   % r1..rp
    a = R \ rvec;                        % AR(p) coefficients

    % One-step prediction errors
    e = zeros(N - p, 1);
    for n = p+1:N
        xhat = a' * xz(n-1:-1:n-p);
        e(n - p) = xz(n) - xhat;
    end

    % Loss function Ep: cumulative squared error (SSE)
    Ep(p) = sum(e.^2);

    MDL(p) = log(Ep(p)) + (p * log(N)) / N;
    AIC(p) = log(Ep(p)) + (2 * p) / N;
    AICc(p) = AIC(p) + (2 * p * (p + 1)) / (N - p - 1);
end

figure;
plot(1:maxp, MDL, '-o', 'DisplayName', 'MDL');
hold on;
plot(1:maxp, AIC, '-s', 'DisplayName', 'AIC');
plot(1:maxp, AICc, '-d', 'DisplayName', 'AICc');
hold off;
xlabel('Model order p', 'FontSize', 14);
ylabel('Criterion value',  'FontSize', 14);
title('Order selection for standardised sunspot data',  'FontSize', 14);
legend('show', 'Location', 'northeast', 'FontSize', 10);
grid on;