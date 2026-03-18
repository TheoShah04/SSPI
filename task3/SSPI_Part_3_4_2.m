clear; clc;

Fs = 32768;  
toneDur = 0.25;       
idleDur = 0.25;   
Nseg = round(toneDur * Fs); % 8192 samples per tone/idle block
Nfft = Nseg; % one FFT per block

rng('shuffle');
randDigits = randi([0, 9], 1, 8);
digitsDial = [0 2 0 randDigits]; % 11 keys
numStr = sprintf('020 %d%d%d%d %d%d%d%d', randDigits(1), randDigits(2), ...
    randDigits(3), randDigits(4), randDigits(5), randDigits(6), ...
    randDigits(7), randDigits(8));
fprintf('Random London number: %s\n', numStr);

y = [];
blockLabel = {};
blockKey = [];
for k = 1:numel(digitsDial)
    d = digitsDial(k);
    [f1, f2] = dtmf_freqs(d);
    n = 0:(Nseg - 1);
    yTone = sin(2*pi*f1*n/Fs) + sin(2*pi*f2*n/Fs);

    y = [y, yTone]; %#ok<AGROW>
    blockLabel{end+1} = sprintf('Key %d', d); %#ok<AGROW>
    blockKey(end+1) = d; %#ok<AGROW>

    if k < numel(digitsDial)
        y = [y, zeros(1, Nseg)]; %#ok<AGROW>
        blockLabel{end+1} = 'Idle'; %#ok<AGROW>
        blockKey(end+1) = NaN; %#ok<AGROW>
    end
end
nBlocks = numel(blockLabel); % should be 21
fprintf('Total blocks: %d (expected 21)\n', nBlocks);

% - Hann window to reduce leakage
win = hann(Nseg, 'periodic');
noverlap = 0;
[S, F, T] = spectrogram(y, win, noverlap, Nfft, Fs);
P = abs(S).^2;
PdB = 10*log10(P + eps);

figure;
imagesc(T, F, PdB);
axis xy;
ylim([0 2000]);
grid on;
colormap(turbo(256));
colorbar;
xlabel('Time (s)', 'FontSize', 14);
ylabel('Frequency (Hz)', 'FontSize', 14);
title('Spectrogram of DTMF Sequence (Hann, Non-overlapping 0.25 s Frames)', 'FontSize', 14);
set(gca, 'FontSize', 14);

figure;
surf(T, F, PdB, 'EdgeColor', 'none');
view(45, 60);
ylim([0 2000]);
zlim([-300 max(PdB(:))]);
colormap(turbo(256));
colorbar;
xlabel('Time (s)', 'FontSize', 14);
ylabel('Frequency (Hz)', 'FontSize', 14);
zlabel('Power (dB)', 'FontSize', 14);
title('3D Spectrogram (Non-Normalised)', 'FontSize', 14);
set(gca, 'FontSize', 14);

Yblk = reshape(y, Nseg, []);
f = (0:(Nfft/2)) * Fs / Nfft;
Pblk = zeros(numel(f), nBlocks);
for b = 1:nBlocks
    Xb = fft(win .* Yblk(:, b), Nfft);
    Pblk(:, b) = abs(Xb(1:(Nfft/2 + 1))).^2;
end

keyIdx = find(~isnan(blockKey));
pick = keyIdx(randperm(numel(keyIdx), 2));
idxA = pick(1);
idxB = pick(2);
idxI = find(isnan(blockKey), 1, 'first');

fprintf('Selected key blocks for PSD plot: %d (%s), %d (%s)\n', ...
    idxA, blockLabel{idxA}, idxB, blockLabel{idxB});

figure;
subplot(3,1,1);
plot(f, 10*log10(Pblk(:, idxA) + eps), 'b', 'LineWidth', 1.1);
grid on; xlim([0 2000]);
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Power (dB)', 'FontSize', 12);
title(sprintf('FFT Segment %d: %s', idxA, blockLabel{idxA}), 'FontSize', 12);
set(gca, 'FontSize', 12);

subplot(3,1,2);
plot(f, 10*log10(Pblk(:, idxI) + eps), 'k', 'LineWidth', 1.1);
grid on; xlim([0 2000]);
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Power (dB)', 'FontSize', 12);
title(sprintf('FFT Segment %d: %s', idxI, blockLabel{idxI}), 'FontSize', 12);
set(gca, 'FontSize', 12);

subplot(3,1,3);
plot(f, 10*log10(Pblk(:, idxB) + eps), 'r', 'LineWidth', 1.1);
grid on; xlim([0 2000]);
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Power (dB)', 'FontSize', 12);
title(sprintf('FFT Segment %d: %s', idxB, blockLabel{idxB}), 'FontSize', 12);
set(gca, 'FontSize', 12);

% Basic analysis printout: two dominant frequencies for tone blocks
fprintf('\nDetected dominant frequencies per tone block:\n');
for b = 1:nBlocks
    if isnan(blockKey(b))
        continue;
    end
    mag = Pblk(:, b);
    mag(1) = 0; % ignore DC
    [~, ord] = sort(mag, 'descend');
    fA = f(ord(1));
    fB = NaN;
    for k = 2:numel(ord)
        if abs(f(ord(k)) - fA) > 40 % enforce separation between two tones
            fB = f(ord(k));
            break;
        end
    end
    [fLo, fHi] = dtmf_freqs(blockKey(b));
    fprintf('Block %2d (%s): peaks ~ %.1f Hz, %.1f Hz | expected %d Hz, %d Hz\n', ...
        b, blockLabel{b}, fA, fB, fLo, fHi);
end

function [fLow, fHigh] = dtmf_freqs(digit)
switch digit
    case 1, fLow = 697; fHigh = 1209;
    case 2, fLow = 697; fHigh = 1336;
    case 3, fLow = 697; fHigh = 1477;
    case 4, fLow = 770; fHigh = 1209;
    case 5, fLow = 770; fHigh = 1336;
    case 6, fLow = 770; fHigh = 1477;
    case 7, fLow = 852; fHigh = 1209;
    case 8, fLow = 852; fHigh = 1336;
    case 9, fLow = 852; fHigh = 1477;
    case 0, fLow = 941; fHigh = 1336;
    otherwise
        error('Unsupported key: %d', digit);
end
end
