clc;
clear;
close all;

syms s

W1 = 1 / (s + 1);
W2 = 2;
W3 = 5 / (s^2 + 3*s + 2);
W4 = 4;

W_final = ((1 + W1)*(1 + W3) - W3*W4) / (1 + W3 + W3*W4*(1 + 2*W2 + W1*W2));
W_final = simplify(W_final)

syms u e1 e2 e3 y1 y2 y3 y4
eqs = [e1 == u - y4, ...
       e2 == e1 + y1 - y4, ...
       e3 == y2 + e1 - y3, ...
       y1 == W1*e1, y2 == W2*e2, y3 == W3*e3, y4 == W4*y3];
sol = solve(eqs, [e1 e2 e3 y1 y2 y3 y4]);
W_check = simplify(sol.e2 / u);
simplify(W_final - W_check)  

[num, den] = numden(W_final);
num = sym2poly(expand(num))
den = sym2poly(expand(den))

