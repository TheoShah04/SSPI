function s = skewness(x)

x = x(:);
x = x(~isnan(x));
n = numel(x);
if n == 0
    s = NaN;
    return;
end

mu = mean(x);
dx = x - mu;
m2 = mean(dx.^2);
m3 = mean(dx.^3);

if m2 == 0
    s = 0;
else
    s = m3 / (m2^(3/2));
end
