
clear;
clc;
close all;

syms s;

W1 = 1 / (s + 1);
W2 = 2;
W3 = 5 / (s^2 + 3*s + 2);
W4 = 4;

W_eq = W3 * (1 + W1*W2) / (1 + W3*W4*(1 + W1*W2) + W3*W2);

[num, den] = numden(W_eq);

num = expand(num);
den = expand(den);

disp('ЭКВИВАЛЕНТНАЯ ПЕРЕДАТОЧНАЯ ФУНКЦИЯ:');

pretty(num / den)

[num, den] = numden(W_eq);

num = double(sym2poly(expand(num)));
den = double(sym2poly(expand(den)));

disp('КОЭФФИЦИЕНТЫ ПЕРЕДАТОЧНОЙ ФУНКЦИИ');
fprintf('Числитель:');
disp(num);
fprintf('Знаменатель:');
disp(den);
