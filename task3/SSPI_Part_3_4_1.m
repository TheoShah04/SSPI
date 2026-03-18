clear; clc;

Fs = 32768;
toneDur = 0.25;
idleDur = 0.25;
N_tone = round(toneDur * Fs);
N_idle = round(idleDur * Fs);

% Generate random number 020 XXXX XXXX
rng('shuffle');
randDigits = randi([0, 9], 1, 8);
numStr = sprintf('020 %d%d%d%d %d%d%d%d', randDigits(1), randDigits(2), ...
    randDigits(3), randDigits(4), randDigits(5), randDigits(6), ...
    randDigits(7), randDigits(8));
digitsDial = [0 2 0 randDigits];

fprintf('Random London number: %s\n', numStr);

% Build full dialed sequence:
y = [];
for k = 1:numel(digitsDial)
    d = digitsDial(k);
    [f1, f2] = dtmf_freqs(d);
    n = 0:(N_tone - 1);
    yTone = sin(2 * pi * f1 * n / Fs) + sin(2 * pi * f2 * n / Fs);
    y = [y, yTone]; %#ok<AGROW>

    if k < numel(digitsDial)
        y = [y, zeros(1, N_idle)]; %#ok<AGROW>
    end
end

T_total = numel(y) / Fs;
fprintf('Dialed sequence duration: %.2f s\n', T_total);
fprintf('Expected duration: %.2f s\n', 0.25 * 21);

% Full sequence plot
t = (0:numel(y)-1) / Fs;
figure;
plot(t, y, 'k');
grid on;
xlabel('Time (s)', 'FontSize', 14);
ylabel('Amplitude', 'FontSize', 14);
title(sprintf('DTMF Sequence for %s', numStr), 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([0 T_total]);

% Plot two different keys and idle sequence (zoomed)
keyA = 2;
keyB = 9;
[f1A, f2A] = dtmf_freqs(keyA);
[f1B, f2B] = dtmf_freqs(keyB);
n = 0:(N_tone - 1);
yA = sin(2 * pi * f1A * n / Fs) + sin(2 * pi * f2A * n / Fs);
yB = sin(2 * pi * f1B * n / Fs) + sin(2 * pi * f2B * n / Fs);
yIdle = zeros(1, N_tone);
tZoom = n / Fs;
tWin = 0.015; % 15ms window
figure;
subplot(3,1,1);
plot(tZoom, yA, 'b');
grid on;
xlabel('Time (s)', 'FontSize', 12);
ylabel('Amp.', 'FontSize', 12);
title(sprintf('Key %d: %d Hz + %d Hz', keyA, f1A, f2A), 'FontSize', 12);
xlim([0 tWin]);
ylim([-2.2 2.2]);
set(gca, 'FontSize', 12);

subplot(3,1,2);
plot(tZoom, yB, 'r');
grid on;
xlabel('Time (s)', 'FontSize', 12);
ylabel('Amp.', 'FontSize', 12);
title(sprintf('Key %d: %d Hz + %d Hz', keyB, f1B, f2B), 'FontSize', 12);
xlim([0 tWin]);
ylim([-2.2 2.2]);
set(gca, 'FontSize', 12);

subplot(3,1,3);
plot(tZoom, yIdle, 'k');
grid on;
xlabel('Time (s)', 'FontSize', 12);
ylabel('Amp.', 'FontSize', 12);
title('Idle sequence (silence)', 'FontSize', 12);
xlim([0 tWin]);
ylim([-2.2 2.2]);
set(gca, 'FontSize', 12);

function [fLow, fHigh] = dtmf_freqs(digit)
    % DTMF frequencies for numeric keys
    switch digit
        case 1
            fLow = 697; fHigh = 1209;
        case 2
            fLow = 697; fHigh = 1336;
        case 3
            fLow = 697; fHigh = 1477;
        case 4
            fLow = 770; fHigh = 1209;
        case 5
            fLow = 770; fHigh = 1336;
        case 6
            fLow = 770; fHigh = 1477;
        case 7
            fLow = 852; fHigh = 1209;
        case 8
            fLow = 852; fHigh = 1336;
        case 9
            fLow = 852; fHigh = 1477;
        case 0
            fLow = 941; fHigh = 1336;
        otherwise
            error('Unsupported key: %d', digit);
    end
end
