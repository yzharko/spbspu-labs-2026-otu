
clc;
clear;

syms s

W1 = 1 / (s + 1);
W2 = 2;
W3 = 5 / (s^2 + 3*s + 2);
W4 = 4;

W_final = W3 * (1 + W1*W2) / (1 + W2*W3 + W3*W4*(1 + W1*W2));

W_final = simplify(W_final);

disp('Эквивалентная передаточная функция:');
disp(W_final);

[num, den] = numden(W_final);

num_coeff = sym2poly(expand(num));
den_coeff = sym2poly(expand(den));

disp('Numerator:');
disp(num_coeff);

disp('Denominator:');
disp(den_coeff);