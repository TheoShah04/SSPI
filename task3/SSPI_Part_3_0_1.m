clear; clc;

rng(13); % reproducible
Ns = [128, 256, 512, 100000];

for i = 1:numel(Ns)
    N = Ns(i);
    x = randn(N, 1);   
    P = pgm(x);
    f = (0:N-1) / N;       

    figure;
    plot(f, P, 'k');
    title(sprintf('Periodogram of WGN (N=%d)', N));
    xlabel('Normalized frequency');
    ylabel('P_b_x(f)');
    grid on;

    fprintf('N=%d: mean(P)=%.3f, var(P)=%.3f\n', N, mean(P), var(P));
end