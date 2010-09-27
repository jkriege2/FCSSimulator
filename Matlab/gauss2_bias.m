function ret=gauss2_bias(b, X)
global sigma1;
ret=b(1)+b(2)*exp(-0.5*((X-b(3)).^2)./(sigma1.^2))+b(4)*exp(-0.5*((X-b(3)).^2)./(b(5).^2));

% b= [offset12 amplitude1 position12 amplitude2 width2]