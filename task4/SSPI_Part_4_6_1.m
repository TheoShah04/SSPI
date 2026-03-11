clear; clc; close all;

mu = 0.05;
N = 4000;
arA = [1, 0.9, 0.2];

rng(67);
eta = randn(N, 1);
x_ar = filter(1, arA, eta);
wTrue = -arA(2:end).';

[W_lms, e_lms] = lms_variant(x_ar, 2, mu, 'lms');
[W_se,  e_se]  = lms_variant(x_ar, 2, mu, 'sign-error');
[W_sr,  e_sr]  = lms_variant(x_ar, 2, mu, 'sign-regressor');
[W_ss,  e_ss]  = lms_variant(x_ar, 2, mu, 'sign-sign');

tail = 200;
fprintf('AR(2) expected coefficients: a1=%.3f, a2=%.3f\n', wTrue(1), wTrue(2));
fprintf('LMS final:        a1=%.4f, a2=%.4f\n', W_lms(1,end), W_lms(2,end));
fprintf('LMS tail mean:    a1=%.4f, a2=%.4f\n', mean(W_lms(1,end-tail+1:end)), mean(W_lms(2,end-tail+1:end)));
fprintf('Sign-error final: a1=%.4f, a2=%.4f\n', W_se(1,end), W_se(2,end));
fprintf('Sign-reg final:   a1=%.4f, a2=%.4f\n', W_sr(1,end), W_sr(2,end));
fprintf('Sign-sign final:  a1=%.4f, a2=%.4f\n\n', W_ss(1,end), W_ss(2,end));

err_lms = norm(W_lms(:,end) - wTrue) / norm(wTrue);
err_se  = norm(W_se(:,end) - wTrue) / norm(wTrue);
err_sr  = norm(W_sr(:,end) - wTrue) / norm(wTrue);
err_ss  = norm(W_ss(:,end) - wTrue) / norm(wTrue);
rt_lms = rise_time_from_error(e_lms, 200, 0.05);
rt_se  = rise_time_from_error(e_se, 200, 0.05);
rt_sr  = rise_time_from_error(e_sr, 200, 0.05);
rt_ss  = rise_time_from_error(e_ss, 200, 0.05);
fprintf('AR(2) relative coefficient error:\n');
fprintf('  LMS          : %.6f\n', err_lms);
fprintf('  sign-error   : %.6f\n', err_se);
fprintf('  sign-regressor: %.6f\n', err_sr);
fprintf('  sign-sign    : %.6f\n\n', err_ss);
fprintf('AR(2) rise time (samples, within 5%% of steady-state MSE):\n');
fprintf('  LMS          : %d\n', rt_lms);
fprintf('  sign-error   : %d\n', rt_se);
fprintf('  sign-regressor: %d\n', rt_sr);
fprintf('  sign-sign    : %d\n\n', rt_ss);

figure('Name', 'AR(2): Coefficient convergence');
tlc = tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile;
plot(ax1, W_lms(1, :), 'LineWidth', 1.4, 'DisplayName', 'LMS');
hold(ax1, 'on');
plot(ax1, W_se(1, :), '--', 'LineWidth', 1.2, 'DisplayName', 'sign-error');
plot(ax1, W_sr(1, :), ':', 'LineWidth', 1.2, 'DisplayName', 'sign-regressor');
plot(ax1, W_ss(1, :), '-.', 'LineWidth', 1.2, 'DisplayName', 'sign-sign');
yline(ax1, wTrue(1), 'k--', 'HandleVisibility', 'off');
hold(ax1, 'off');
grid(ax1, 'on');
xlabel(ax1, 'n', 'FontSize', 14);
ylabel(ax1, 'a_1[n]', 'FontSize', 14);
title(ax1, 'Coefficient a_1 convergence', 'FontSize', 14);
legend(ax1, 'Location', 'best');

ax2 = nexttile;
plot(ax2, W_lms(2, :), 'LineWidth', 1.4, 'DisplayName', 'LMS');
hold(ax2, 'on');
plot(ax2, W_se(2, :), '--', 'LineWidth', 1.2, 'DisplayName', 'sign-error');
plot(ax2, W_sr(2, :), ':', 'LineWidth', 1.2, 'DisplayName', 'sign-regressor');
plot(ax2, W_ss(2, :), '-.', 'LineWidth', 1.2, 'DisplayName', 'sign-sign');
yline(ax2, wTrue(2), 'k--', 'HandleVisibility', 'off');
hold(ax2, 'off');
grid(ax2, 'on');
xlabel(ax2, 'n', 'FontSize', 14);
ylabel(ax2, 'a_2[n]', 'FontSize', 14);
title(ax2, 'Coefficient a_2 convergence', 'FontSize', 14);
legend(ax2, 'Location', 'best');

sgtitle(tlc, 'AR(2) identification: LMS vs sign LMS variants', 'FontSize', 14);

% Error curves (AR case)
figure('Name', 'AR(2): Learning curves');
plot(movmean(e_lms.^2, 50), 'LineWidth', 1.2, 'DisplayName', 'LMS');
hold on;
plot(movmean(e_se.^2, 50), 'LineWidth', 1.2, 'DisplayName', 'sign-error');
plot(movmean(e_sr.^2, 50), 'LineWidth', 1.2, 'DisplayName', 'sign-regressor');
plot(movmean(e_ss.^2, 50), 'LineWidth', 1.2, 'DisplayName', 'sign-sign');
hold off;
grid on;
xlabel('n', 'FontSize', 14);
ylabel('Smoothed e^2[n]', 'FontSize', 14);
title('AR(2) learning curves', 'FontSize', 14);
legend('Location', 'best');

fs = 44100;
Nsp = 1000;
audioDir = fullfile(fileparts(mfilename('fullpath')), 'audio');
if ~exist(audioDir, 'dir')
    audioDir = fileparts(mfilename('fullpath'));
end
speechFile = fullfile(audioDir, 'e.wav');

[x_sp, ok] = load_audio_segment(speechFile, fs, Nsp);
if ok
    p = 16;  % predictor order
    [Wl, el] = lms_variant(x_sp, p, mu, 'lms');
    [Wse, ese] = lms_variant(x_sp, p, mu, 'sign-error');
    [Wsr, esr] = lms_variant(x_sp, p, mu, 'sign-regressor');
    [Wss, ess] = lms_variant(x_sp, p, mu, 'sign-sign');

    % Print error of estimate (speech) as steady-state MSE
    win = 200;
    [~, speechName, speechExt] = fileparts(speechFile);
    fprintf('Speech "%s%s" steady-state MSE (last %d samples):\n', speechName, speechExt, win);
    fprintf('  LMS          : %.6f\n', mean(el(end-win+1:end).^2));
    fprintf('  sign-error   : %.6f\n', mean(ese(end-win+1:end).^2));
    fprintf('  sign-regressor: %.6f\n', mean(esr(end-win+1:end).^2));
    fprintf('  sign-sign    : %.6f\n\n', mean(ess(end-win+1:end).^2));

    rt_lms_sp = rise_time_from_error(el, win, 0.05);
    rt_se_sp  = rise_time_from_error(ese, win, 0.05);
    rt_sr_sp  = rise_time_from_error(esr, win, 0.05);
    rt_ss_sp  = rise_time_from_error(ess, win, 0.05);
    fprintf('Speech "%s%s" rise time (samples, within 5%% of steady-state MSE):\n', speechName, speechExt);
    fprintf('  LMS          : %d\n', rt_lms_sp);
    fprintf('  sign-error   : %d\n', rt_se_sp);
    fprintf('  sign-regressor: %d\n', rt_sr_sp);
    fprintf('  sign-sign    : %d\n\n', rt_ss_sp);

    % Learning curves (speech)
    figure('Name', 'Speech: Learning curves');
    plot(movmean(el.^2, 50), 'LineWidth', 1.2, 'DisplayName', 'LMS');
    hold on;
    plot(movmean(ese.^2, 50), 'LineWidth', 1.2, 'DisplayName', 'sign-error');
    plot(movmean(esr.^2, 50), 'LineWidth', 1.2, 'DisplayName', 'sign-regressor');
    plot(movmean(ess.^2, 50), 'LineWidth', 1.2, 'DisplayName', 'sign-sign');
    hold off;
    grid on;
    xlabel('n', 'FontSize', 14);
    ylabel('Smoothed e^2[n]', 'FontSize', 14);
    title('Speech learning curves ("a.wav")', 'FontSize', 14);
    legend('Location', 'best');

    % Convergence of first two coefficients (speech)
    figure('Name', 'Speech: Coefficient convergence (first two taps)');
    plot(Wl(1, :), 'LineWidth', 1.1, 'DisplayName', 'LMS w_1');
    hold on;
    plot(Wl(2, :), 'LineWidth', 1.1, 'DisplayName', 'LMS w_2');
    plot(Wse(1, :), '--', 'LineWidth', 1.1, 'DisplayName', 'sign-error w_1');
    plot(Wse(2, :), '--', 'LineWidth', 1.1, 'DisplayName', 'sign-error w_2');
    plot(Wsr(1, :), ':', 'LineWidth', 1.1, 'DisplayName', 'sign-regressor w_1');
    plot(Wsr(2, :), ':', 'LineWidth', 1.1, 'DisplayName', 'sign-regressor w_2');
    plot(Wss(1, :), '-.', 'LineWidth', 1.1, 'DisplayName', 'sign-sign w_1');
    plot(Wss(2, :), '-.', 'LineWidth', 1.1, 'DisplayName', 'sign-sign w_2');
    hold off;
    grid on;
    xlabel('n', 'FontSize', 14);
    ylabel('Coefficient value', 'FontSize', 14);
    title('Speech coefficients (first two taps)', 'FontSize', 14);
    legend('Location', 'best');
else
    fprintf('Speech file not found: %s\n', speechFile);
end

function [W, e] = lms_variant(x, p, mu, mode)
    N = numel(x);
    W = zeros(p, N);
    e = zeros(N, 1);
    w = zeros(p, 1);
    for n = p+1:N
        xvec = x(n-1:-1:n-p);
        yhat = w.' * xvec;
        e(n) = x(n) - yhat;
        switch mode
            case 'lms'
                w = w + mu * e(n) * xvec;
            case 'sign-error'
                w = w + mu * sign(e(n)) * xvec;
            case 'sign-regressor'
                w = w + mu * e(n) * sign(xvec);
            case 'sign-sign'
                w = w + mu * sign(e(n)) .* sign(xvec);
        end
        W(:, n) = w;
    end
end

function [x, ok] = load_audio_segment(filePath, fs, N)
    ok = false;
    if ~exist(filePath, 'file')
        x = [];
        return;
    end
    [x, fsIn] = audioread(filePath);
    if size(x, 2) > 1
        x = mean(x, 2);
    end
    if fsIn ~= fs
        if exist('resample', 'file') == 2
            x = resample(x, fs, fsIn);
        else
            tOld = (0:numel(x)-1) / fsIn;
            tNew = (0:round(numel(x) * fs / fsIn) - 1) / fs;
            x = interp1(tOld, x, tNew, 'linear', 'extrap');
        end
    end
    if numel(x) < N
        x = [x; zeros(N - numel(x), 1)];
    elseif numel(x) > N
        win = N;
        energy = conv(x.^2, ones(win, 1), 'valid');
        [~, idx] = max(energy);
        x = x(idx:idx+N-1);
    end
    x = x - mean(x);
    ok = true;
end

function rt = rise_time_from_error(e, win, tolFrac)
    % Rise time: first index where smoothed MSE is within (1+tolFrac) of steady-state MSE
    if nargin < 2, win = 200; end
    if nargin < 3, tolFrac = 0.10; end
    e2 = e.^2;
    mse_ss = mean(e2(end-win+1:end));
    thr = (1 + tolFrac) * mse_ss;
    mse_sm = movmean(e2, 50);
    idx = find(mse_sm <= thr, 1, 'first');
    if isempty(idx)
        rt = numel(e);
    else
        rt = idx;
    end
end
