function ret=gauss_bias(b, X)

ret=b(1)+b(2)*exp(-0.5*((X-b(3)).^2)./(b(4).^2));