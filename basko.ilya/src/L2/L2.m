clc; clear; close all;

syms s;

W1 = 1/(s + 1);
W2 = 2;
W3 = 5/(s^2 + 3*s + 2);
W4 = 4;

A = 1 + W1*W2;
W_inner = W3/(1 + W3*W2*W4);
W_eq = simplifyFraction(A*W_inner/(1 + A*W_inner*W4), 'Expand', true)
