clear; clc;

Fs = 32768;
toneDur = 0.25;
idleDur = 0.25;
Nseg = round(toneDur * Fs);
Nfft = Nseg;
win = hann(Nseg, 'periodic');
noverlap = 0;

rng('shuffle');
randDigits = randi([0, 9], 1, 8);
digitsDial = [0 2 0 randDigits];
numStr = sprintf('020 %d%d%d%d %d%d%d%d', randDigits(1), randDigits(2), ...
    randDigits(3), randDigits(4), randDigits(5), randDigits(6), ...
    randDigits(7), randDigits(8));
[yClean, blockLabel, blockKey] = build_dtmf_sequence(digitsDial, Fs, Nseg);
nBlocks = numel(blockLabel);
sigVar = var(yClean);

noiseScale = [0.01, 0.25, 4.00];
noiseVar = sigVar * noiseScale;
snrDB = 10 * log10(sigVar ./ noiseVar);

fprintf('Random London number: %s\n', numStr);
for c = 1:numel(noiseVar)
    fprintf('Case %d: noise variance = %.4f, approx SNR = %.2f dB\n', ...
        c, noiseVar(c), snrDB(c));
end

% 2D spectrograms for the three noisy cases
figure;
for c = 1:numel(noiseVar)
    yNoisy = yClean + sqrt(noiseVar(c)) * randn(size(yClean));
    [S, F, T] = spectrogram(yNoisy, win, noverlap, Nfft, Fs);
    PdB = 10*log10(abs(S).^2 + eps);

    subplot(3,1,c);
    imagesc(T, F, PdB);
    axis xy;
    ylim([0 2000]);
    colormap(turbo(256));
    colorbar;
    grid on;
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('Frequency (Hz)', 'FontSize', 12);
    title(sprintf('Noisy Spectrogram Case %d (\\sigma_v^2=%.3f, SNR=%.2f dB)', ...
        c, noiseVar(c), snrDB(c)), 'FontSize', 12);
    set(gca, 'FontSize', 12);
end

% 3D non-normalized spectrogram for the noisiest case
yNoisyWorst = yClean + sqrt(noiseVar(end)) * randn(size(yClean));
[S3, F3, T3] = spectrogram(yNoisyWorst, win, noverlap, Nfft, Fs);
PdB3 = 10*log10(abs(S3).^2 + eps);
PdB3_pos = max(PdB3, 20);
figure;
surf(T3, F3, PdB3_pos, 'EdgeColor', 'none');
view(45, 60);
ylim([0 2000]);
zlim([0 max(PdB3_pos(:))]);
colormap(turbo(256));
colorbar;
grid on;
xlabel('Time (s)', 'FontSize', 14);
ylabel('Frequency (Hz)', 'FontSize', 14);
zlabel('Power (dB)', 'FontSize', 14);
title('3D Spectrogram (Noisiest Case, Non-Normalized)', 'FontSize', 14);
set(gca, 'FontSize', 14);

% FFT segments for two random keys + one idle, for each noise case
keyIdx = find(~isnan(blockKey));
pick = keyIdx(randperm(numel(keyIdx), 2));
idleIdx = find(isnan(blockKey), 1, 'first');
showIdx = [pick(1), idleIdx, pick(2)];
f = (0:(Nfft/2)) * Fs / Nfft;

for c = 1:numel(noiseVar)
    yNoisy = yClean + sqrt(noiseVar(c)) * randn(size(yClean));
    Yblk = reshape(yNoisy, Nseg, []);
    Pblk = zeros(numel(f), nBlocks);
    for b = 1:nBlocks
        Xb = fft(win .* Yblk(:, b), Nfft);
        Pblk(:, b) = abs(Xb(1:(Nfft/2 + 1))).^2;
    end

    figure;
    for r = 1:3
        subplot(3,1,r);
        plot(f, 10*log10(Pblk(:, showIdx(r)) + eps), 'LineWidth', 1.1);
        grid on;
        xlim([0 2000]);
        xlabel('Frequency (Hz)', 'FontSize', 12);
        ylabel('Power (dB)', 'FontSize', 12);
        title(sprintf('Case %d FFT Segment %d: %s', c, showIdx(r), blockLabel{showIdx(r)}), ...
            'FontSize', 12);
        set(gca, 'FontSize', 12);
    end
end

% Single-run tone identification on this generated sequence
fprintf('\nSingle-run tone identification accuracy:\n');
for c = 1:numel(noiseVar)
    yNoisy = yClean + sqrt(noiseVar(c)) * randn(size(yClean));
    acc = tone_ident_accuracy(yNoisy, blockKey, Fs, Nseg, win, Nfft);
    fprintf('Case %d: accuracy = %.2f %%\n', c, 100*acc);
end

nTrials = 60;
accMC = zeros(numel(noiseVar), nTrials);
for tr = 1:nTrials
    rd = randi([0, 9], 1, 8);
    dig = [0 2 0 rd];
    [y0, ~, keys0] = build_dtmf_sequence(dig, Fs, Nseg);
    sigVar0 = var(y0);
    nv0 = sigVar0 * noiseScale;
    for c = 1:numel(nv0)
        yn = y0 + sqrt(nv0(c)) * randn(size(y0));
        accMC(c, tr) = tone_ident_accuracy(yn, keys0, Fs, Nseg, win, Nfft);
    end
end

accMean = mean(accMC, 2);
accStd = std(accMC, 0, 2);

fprintf('\nMonte Carlo (%d trials) tone-identification results:\n', nTrials);
for c = 1:numel(noiseVar)
    fprintf('Case %d (SNR %.2f dB): mean accuracy = %.2f %% (std %.2f %%)\n', ...
        c, snrDB(c), 100*accMean(c), 100*accStd(c));
end

function acc = tone_ident_accuracy(y, blockKey, Fs, Nseg, win, Nfft)
    Yblk = reshape(y, Nseg, []);
    toneIdx = find(~isnan(blockKey));
    nCorrect = 0;
    tolHz = 20;
    for m = 1:numel(toneIdx)
        b = toneIdx(m);
        d = blockKey(b);
        xb = win .* Yblk(:, b);
        X = fft(xb, Nfft);
        f = (0:(Nfft/2)) * Fs / Nfft;
        P = abs(X(1:(Nfft/2 + 1))).^2;
        dB = 10*log10(P + eps);
        band = (f >= 600 & f <= 1600);
        fBand = f(band);
        dBand = dB(band);

        prom = max(3, 0.15 * (max(dBand) - median(dBand)));
        [pk, loc] = findpeaks(dBand, fBand, ...
            'MinPeakProminence', prom, ...
            'MinPeakDistance', 60);

        if numel(pk) < 2
            continue;
        end
        [~, ord] = sort(pk, 'descend');
        fDet = sort([loc(ord(1)), loc(ord(2))]);
        [fLo, fHi] = dtmf_freqs(d);
        fExp = sort([fLo, fHi]);

        if all(abs(fDet - fExp) <= tolHz)
            nCorrect = nCorrect + 1;
        end
    end
    acc = nCorrect / numel(toneIdx);
end

function [y, blockLabel, blockKey] = build_dtmf_sequence(digitsDial, Fs, Nseg)
    y = [];
    blockLabel = {};
    blockKey = [];
    for k = 1:numel(digitsDial)
        d = digitsDial(k);
        [f1, f2] = dtmf_freqs(d);
        n = 0:(Nseg-1);
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
