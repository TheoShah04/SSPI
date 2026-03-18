clear; clc;

Fs = 32768;                
toneDur = 0.25;          
idleDur = 0.25;          
N_tone = round(toneDur * Fs);
N_idle = round(idleDur * Fs);

rng('shuffle');
randDigits = randi([0, 9], 1, 8);
digitsDial = [0 2 0 randDigits];
numStr = sprintf('020 %d%d%d%d %d%d%d%d', randDigits(1), randDigits(2), ...
    randDigits(3), randDigits(4), randDigits(5), randDigits(6), ...
    randDigits(7), randDigits(8));

[yClean, blockLabel, blockKey] = build_dtmf_sequence(digitsDial, Fs, N_tone, N_idle);
t = (0:numel(yClean)-1) / Fs;
sigVar = var(yClean);

% Three noise conditions
noiseScale = [0.01, 0.25, 4.00];
noiseVar = sigVar * noiseScale;
snrDB = 10 * log10(sigVar ./ noiseVar);

fprintf('Random London number: %s\n', numStr);
fprintf('Sequence duration: %.2f s (expected %.2f s)\n', numel(yClean)/Fs, 0.25*21);
for k = 1:numel(noiseVar)
    fprintf('Case %d: noise variance = %.4f, approx SNR = %.2f dB\n', ...
        k, noiseVar(k), snrDB(k));
end

yNoisy = cell(1, numel(noiseVar));
for k = 1:numel(noiseVar)
    yNoisy{k} = yClean + sqrt(noiseVar(k)) * randn(size(yClean));
end

figure;
subplot(4,1,1);
plot(t, yClean, 'k');
grid on;
xlim([0 t(end)]);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Amp.', 'FontSize', 12);
title(sprintf('Clean DTMF Sequence: %s', numStr), 'FontSize', 12);
set(gca, 'FontSize', 12);

for k = 1:numel(noiseVar)
    subplot(4,1,k+1);
    plot(t, yNoisy{k}, 'LineWidth', 0.8);
    grid on;
    xlim([0 t(end)]);
    xlabel('Time (s)', 'FontSize', 12);
    ylabel('Amp.', 'FontSize', 12);
    title(sprintf('Noisy Sequence Case %d: \\sigma_v^2=%.3f, SNR=%.2f dB', ...
        k, noiseVar(k), snrDB(k)), 'FontSize', 12);
    set(gca, 'FontSize', 12);
end

keyIdx = find(~isnan(blockKey));
pick = keyIdx(randperm(numel(keyIdx), 2));
idleIdx = find(isnan(blockKey), 1, 'first');
showIdx = [pick(1), idleIdx, pick(2)];
showNames = {blockLabel{showIdx(1)}, blockLabel{showIdx(2)}, blockLabel{showIdx(3)}};

YblkClean = reshape(yClean, N_tone, []);
YblkNoisy = cellfun(@(yy) reshape(yy, N_tone, []), yNoisy, 'UniformOutput', false);
n = (0:N_tone-1) / Fs;
tWin = 0.015; % 15ms

figure;
for r = 1:3
    for c = 1:4
        subplot(3,4,(r-1)*4 + c);
        if c == 1
            plot(n, YblkClean(:, showIdx(r)), 'k', 'LineWidth', 0.9);
            ttl = sprintf('%s (Clean)', showNames{r});
        else
            plot(n, YblkNoisy{c-1}(:, showIdx(r)), 'LineWidth', 0.9);
            ttl = sprintf('%s (Case %d)', showNames{r}, c-1);
        end
        grid on;
        xlim([0 tWin]);
        ylim([-4 4]);
        xlabel('Time (s)', 'FontSize', 10);
        ylabel('Amp.', 'FontSize', 10);
        title(ttl, 'FontSize', 10);
        set(gca, 'FontSize', 10);
    end
end

function [y, blockLabel, blockKey] = build_dtmf_sequence(digitsDial, Fs, N_tone, N_idle)
    y = [];
    blockLabel = {};
    blockKey = [];
    for k = 1:numel(digitsDial)
        d = digitsDial(k);
        [f1, f2] = dtmf_freqs(d);
        n = 0:(N_tone-1);
        yTone = sin(2*pi*f1*n/Fs) + sin(2*pi*f2*n/Fs);
        y = [y, yTone]; %#ok<AGROW>
        blockLabel{end+1} = sprintf('Key %d', d); %#ok<AGROW>
        blockKey(end+1) = d; %#ok<AGROW>
        if k < numel(digitsDial)
            y = [y, zeros(1, N_idle)]; %#ok<AGROW>
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
