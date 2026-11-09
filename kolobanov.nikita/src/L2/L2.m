clear;
close all;
clc;

syms s;

W1 = 1 / (s + 1);
W2 = 2;
W3 = 5 / (s^2 + 3 * s + 2);
W4 = 4;

A = 1 + W1 * W2;
B = W3 / (1 + W2 * W3 * W4);
W_final = simplify(A * B / (1 + A * B * W4));

[num, den] = numden(W_final);
numerator = sym2poly(expand(num));
denominator = sym2poly(expand(den));

disp('W(s):');
disp(W_final);
disp('Числитель:');
disp(numerator);
disp('Знаменатель:');
disp(denominator);
