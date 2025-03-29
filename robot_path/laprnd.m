function y  = laprnd(m, n, mu, sigma) % Generisanje Laplasovog šuma

u = rand(m,n) - 0.5;
b = sigma/sqrt(2);
y = mu - b*sign(u).*log(1 - 2*abs(u));
