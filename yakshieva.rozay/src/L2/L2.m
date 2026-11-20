clear;
clc;
close all;

syms s

W1 = 1 / (s + 1);
W2 = 2;
W3 = 5 / (s^2 + 3 * s + 2);
W4 = 4;

W_final = simplify( ...
    (1 + W3 + W2 * W3 * W4) / ...
    (1 + W3 + W2 * W3 * W4 + W3 * W4 * (1 + W2 + W1 * W2)));

W_final

num = [1 4 50 47];
den = [1 4 110 147];
