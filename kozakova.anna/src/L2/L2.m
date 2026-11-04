clc;
clear;
close all;

syms s

W1 = 1/(s + 1);
W2 = 2;
W3 = 5/(s^2 + 3*s + 2);

W = simplify((1 + W3 + W2*W3)/(1 + 2*W3 + W2*W3 + W1*W2*W3))
[num, den] = numden(W);
Numerator = sym2poly(expand(num))
Denominator = sym2poly(expand(den))
