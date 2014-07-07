function F = myfun_lsqnl_070614(x)
k = 1:10;
F = (2 + 2*k-exp(k*x(1))-exp(k*x(2))).^2;


