pkg load symbolic
syms s

W1 = 1 / (s + 1);
W2 = 2;
W3 = 5 / (s^2 + 3*s + 2);
W4 = 4;

W_final = (1 + W3*W4*W2) / (1+W3*W4*(1 + W2 + W1*W2));
simplify(W_final)

