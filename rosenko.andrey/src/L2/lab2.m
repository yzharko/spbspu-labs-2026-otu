syms s w1 w2 w3 w4

B = 1 + w3 + w2*w3*w4;
W_general = B / (B + w3*w4*(1 + w2 + w1*w2));

W1 = 1/(s + 1);
W2 = 2;
W3 = 5/(s^2 + 3*s + 2);
W4 = 4;

W_sub = subs(W_general, [w1 w2 w3 w4], [W1 W2 W3 W4]);
[N, D] = numden(simplify(W_sub));
W = expand(N)/expand(D)

num = sym2poly(N)
den = sym2poly(D)
