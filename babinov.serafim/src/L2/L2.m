clc;
clear;
close all;

syms s

W1 = 1 / (s + 1);
W2 = 2;
W3 = 5 / (s^2 + 3*s + 2);

W_result = (W1 + W1*W3 - W3) / (1 + 2*W3 + W2*W3*(W1 + 1));
W_result = simplify(W_result);
