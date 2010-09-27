function ret=gauss(b, X)

ret=0+b(3)*exp(-0.5*((X-b(1)).^2)./(b(2).^2));