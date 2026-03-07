clear; clc;

rng(42); 
N = 1064;
x = randn(N, 1);

b = 1;
a = [1 0.9];
y_full = filter(b, a, x);

trim = 40;
y = y_full(trim+1:end);
x_trim = x(trim+1:end);
n = (1:numel(x_trim)).';

figure;
subplot(2,1,1);
plot(n, x_trim, 'k');
title('Input WGN x[n] (trimmed)');
xlabel('n');
ylabel('Amplitude');
xlim([1 numel(x_trim)]);
grid on;

subplot(2,1,2);
plot(n, y, 'b');
title('Filtered output y[n] (AR(1), transient removed)');
xlabel('n');
ylabel('Amplitude');
xlim([1 numel(y)]);
grid on;

[h, w] = freqz(b, a, 512);
H2 = abs(h).^2;
target = max(H2) / 2; % -3dB
idx = find(H2 >= target, 1, 'first');
fc = w(idx) / (2*pi);
fprintf('Cutoff frequency (-3 dB) ~ %.4f (normalized)\n', fc);
Py = pgm(y);
Np = numel(Py);
nhalf = floor(Np/2);
fy = (0:nhalf) / Np;

figure;
subplot(2,1,1);
plot(w/(2*pi), abs(h).^2, 'LineWidth', 1.2, 'DisplayName', 'Theoretical PSD');
hold on;
plot(fy, Py(1:nhalf+1), 'r', 'DisplayName', 'Periodogram of y');
hold off;
title('Theoretical PSD of AR(1) Process (Full)');
xlabel('Normalized frequency');
ylabel('$|H(e^{j2\pi f})|^2$', 'Interpreter', 'latex');
legend('show', 'Location', 'northeast');
xlim([0 0.5]);
grid on;

subplot(2,1,2);
plot(w/(2*pi), abs(h).^2, 'LineWidth', 1.2, 'DisplayName', 'Theoretical PSD');
hold on;
plot(fy, Py(1:nhalf+1), 'r', 'DisplayName', 'Periodogram of y');
hold off;
title('Theoretical PSD of AR(1) Process (Zoom)');
xlabel('Normalized frequency');
ylabel('$|H(e^{j2\pi f})|^2$', 'Interpreter', 'latex');
legend('show', 'Location', 'northeast');
xlim([0.4 0.5]);
grid on;

% Model-based PSD estimate using xcorr
r = xcorr(y, 1, 'biased');
R0 = r(2);
R1 = r(3);
a1_hat = -R1 / R0;
sigma2_hat = R0 + a1_hat * R1;

fprintf('Estimated a1 = %.4f, estimated sigma_x^2 = %.4f\n', a1_hat, sigma2_hat);

[h_mb, w_mb] = freqz(1, [1 a1_hat], 512);
P_mb = sigma2_hat * abs(h_mb).^2;

figure;
subplot(2,1,1);
plot(w/(2*pi), H2, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Theoretical PSD');
hold on;
plot(w_mb/(2*pi), P_mb, 'b', 'LineWidth', 1.2, 'DisplayName', 'Model-based PSD');
plot(fy, Py(1:nhalf+1), 'r', 'DisplayName', 'Periodogram');
hold off;
title('Model-based PSD vs Periodogram (Full)', 'FontSize', 16);
xlabel('Normalized frequency', 'FontSize', 14);
ylabel('PSD', 'FontSize', 14);
lgd1 = legend('show', 'Location', 'northeast');
set(lgd1, 'FontSize', 16);
set(gca, 'FontSize', 14);
xlim([0 0.5]);
grid on;

subplot(2,1,2);
plot(w/(2*pi), H2, 'k--', 'LineWidth', 1.2, 'DisplayName', 'Theoretical PSD');
hold on;
plot(w_mb/(2*pi), P_mb, 'b', 'LineWidth', 1.2, 'DisplayName', 'Model-based PSD');
plot(fy, Py(1:nhalf+1), 'r', 'DisplayName', 'Periodogram');
hold off;
title('Model-based PSD vs Periodogram (Zoom)', 'FontSize', 16);
xlabel('Normalized frequency', 'FontSize', 14);
ylabel('PSD', 'FontSize', 14);
lgd2 = legend('show', 'Location', 'northeast');
set(lgd2, 'FontSize', 16);
set(gca, 'FontSize', 14);
xlim([0.4 0.5]);
grid on;
