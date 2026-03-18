clear; clc; close all;

orders = [1, 2, 10, 100];   
Nfft = 1024;

S = load('sunspot.dat');
if isstruct(S)
    fn = fieldnames(S);
    D = S.(fn{1});
else
    D = S;
end

if size(D, 2) >= 2
    year = D(:, 1);
    sunspot = D(:, 2);
else
    year = (1:numel(D)).';
    sunspot = D(:);
end

sunspot = sunspot(:);
year = year(:);

series = {sunspot, sunspot - mean(sunspot)};
seriesLabel = {'Original series', 'Mean-centered series'};

for s = 1:numel(series)
    x = series{s};
    N = numel(x);

    Pp = pgm(x);
    f = (0:N-1) / N;
    idx = 1:(floor(N/2) + 1);
    if s == 2
        Pp(1) = NaN; % remove DC component for mean-centered plot
    end

    figure;
    semilogy(f(idx), Pp(idx), 'Color', [0.6 0.6 0.6], 'DisplayName', 'Periodogram');
    hold on;
    colors = lines(numel(orders));

    for k = 1:numel(orders)
        p = orders(k);
        [aHat, sigma2Hat] = aryule(x, p);
        [h, w] = freqz(sqrt(sigma2Hat), aHat, Nfft);
        Pm = abs(h).^2;
        if s == 2
            Pm(1) = NaN; % remove DC component for mean-centered plot
        end

        semilogy(w/(2*pi), Pm, 'LineWidth', 1.3, 'Color', colors(k, :), ...
            'DisplayName', sprintf('Model-based AR(%d)', p));
    end

    hold off;
    grid on;
    xlim([0 0.5]);
    xlabel('Normalised frequency', 'FontSize', 14);
    ylabel('PSD', 'FontSize', 14);
    title(sprintf('Sunspot PSD: %s', seriesLabel{s}), 'FontSize', 16);
    set(gca, 'FontSize', 14);
    legend('show', 'Location', 'northeast', 'FontSize', 10);

    X = fft(x, N);
    Pfft = (1 / N) * (abs(X).^2);
    f_fft = (0:N-1) / N;
    idx_fft = 1:(floor(N/2) + 1);
    if s == 2
        Pfft(1) = NaN;
    end

    figure;
    semilogy(f_fft(idx_fft), Pfft(idx_fft), 'k', 'LineWidth', 1.1);
    grid on;
    xlim([0 0.5]);
    xlabel('Normalised frequency', 'FontSize', 14);
    ylabel('Periodogram (|FFT|^2 / N)', 'FontSize', 14);
    title(sprintf('FFT-Based Periodogram: %s', seriesLabel{s}), 'FontSize', 16);
    set(gca, 'FontSize', 14);
end

figure;
subplot(2,1,1);
plot(year, sunspot, 'k');
grid on;
xlabel('Year', 'FontSize', 14);
ylabel('Sunspot number', 'FontSize', 14);
title('Sunspot series (original)', 'FontSize', 14);
set(gca, 'FontSize', 14);
ylim_orig = ylim;
ylim_span = diff(ylim_orig);

subplot(2,1,2);
plot(year, sunspot - mean(sunspot), 'b');
grid on;
xlabel('Year', 'FontSize', 14);
ylabel('Sunspot number', 'FontSize', 14);
title('Sunspot series (mean-centered)', 'FontSize', 14);
set(gca, 'FontSize', 14);
ylim([-ylim_span/2, ylim_span/2]); 
