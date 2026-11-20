syms s

W1 = 1/(s + 1);
W2 = 2;
W3 = 5/(s^2 + 3*s + 2);
W4 = 4;

W_final = simplify((W1 - W3*W4) / ...
    (1 + W3*W4 + W2*W3*W4*(1 + W1)));

disp(W_final);
