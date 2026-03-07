function P = pgm(x)

isRow = isrow(x);
x = x(:);
N = numel(x);

X = fft(x, N);
P = (1 / N) * (abs(X).^2);

if isRow
    P = P.';
end
end
