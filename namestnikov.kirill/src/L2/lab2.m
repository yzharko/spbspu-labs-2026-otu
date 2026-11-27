clc;
clear;
close all;

syms s
W1 = 1/(s + 1);
W2 = 2;
W3 = 5/(s^2 + 3*s + 2);
W4 = 4;

G = W3*(1 + W1*W2)/(1 + W2*W3);

W_final = simplify((W1 - G)/(1 + W4*G))
