clc;
clear;
close all;

syms s

% Заданные передаточные функции.
W1 = 1 / (s + 1);
W2 = 2;
W3 = 5 / (s^2 + 3*s + 2);
W4 = 4;

% выход y находится перед W2.
W_final = simplify((W1 - W3*W4) / ...
    (1 + W3*W4*(1 + W2 + W1*W2)))

% Коэффициенты для блока Transfer Fcn.
[num, den] = numden(W_final);
num = sym2poly(expand(num))
den = sym2poly(expand(den))
