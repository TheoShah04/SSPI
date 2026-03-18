N = 1000;
numPairs = 100;
rng(42);

A1 = -2.5 + 5.0 * rand(numPairs, 1);
A2 = -1.5 + 3.0 * rand(numPairs, 1);

stable = false(numPairs, 1);

for i = 1:numPairs
    a1 = A1(i);
    a2 = A2(i);

    % Generate AR(2) process x[n] = a1 x[n-1] + a2 x[n-2] + w[n]
    w = randn(N, 1);
    x = zeros(N, 1);
    for n = 3:N
        x(n) = a1 * x(n-1) + a2 * x(n-2) + w(n);
    end

    r = roots([1 -a1 -a2]);
    stable(i) = all(abs(r) < 1);
end

figure;
plot(A1(stable), A2(stable), 'k*');
hold on;

a1_line = linspace(-2.5, 2.5, 200);
plot(a1_line, 1 - a1_line, 'r--');   
plot(a1_line, 1 + a1_line, 'r--');
plot(a1_line, -1 * ones(size(a1_line)), 'r--');

hold off;
xlabel('a_1', 'FontSize', 14);
ylabel('a_2', 'FontSize', 14);
title('Convergence region for AR(2) coefficients', 'FontSize', 14);
axis([-2.5 2.5 -1.5 1.5]);
grid on;
legend('Convergent (simulated)', 'Stability boundary', 'Location', 'northeast');