clear; clc;

load sunspot.dat
x = sunspot(:, 2);

Ns = [5 20 250];
figure;
h1 = [];
h2 = [];
ax1 = [];
for i = 1:numel(Ns)
    N = Ns(i);
    xi = x(1:N);
    xzi = xi - mean(xi);

    [r, lags] = xcorr(xi, 'unbiased');
    [rz, lagsz] = xcorr(xzi, 'unbiased');

    subplot(3,1,i);
    h1i = stem(lags, r, 'filled', 'DisplayName', 'Original');
    hold on;
    h2i = stem(lagsz, rz, 'filled', 'DisplayName', 'Zero-mean');
    hold off;
    if i == 3
        set([h1i h2i], 'LineWidth', 0.5, 'MarkerSize', 3);
    end
    if i == 1
        h1 = h1i;
        h2 = h2i;
        ax1 = gca;
    end
    xlabel('Lag \tau', 'FontSize', 14);
    ylabel('ACF', 'FontSize', 14);
    title(sprintf('ACF (N = %d)', N), 'FontSize', 14);
    grid on;
end
lgd = legend(ax1, [h1 h2], {'Original','Zero-mean'}, 'Location', 'northeastoutside');
legend boxoff;
set(lgd, 'Units', 'normalized');
pos = get(lgd, 'Position');
pos(1) = 0.70; 
pos(2) = 0.93; 
set(lgd, 'Position', pos);
