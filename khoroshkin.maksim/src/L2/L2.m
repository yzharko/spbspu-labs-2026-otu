clc;
clear;
close all;

syms s

W1 = 1 / (s + 1);
W2 = 2;
W3 = 5 / (s^2 + 3*s + 2);
W4 = 4;

A = 1 + W2 + W1*W2;

W_final = (W3*A) / (1 + W3 + W2*W3*W4 + W3*W4*A);

pretty(simplify(W_final))
