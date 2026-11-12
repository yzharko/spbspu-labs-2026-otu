syms s

W1 = 1/(s + 1);
W2 = 2;
W3 = 5/(s^2 + 3*s + 2);
W4 = 4;

A = 1 + W2 + W1*W2;

W_final = simplify(W3*A / ...
    (1 + W3 + W2*W3*W4 + W3*W4*A));

disp('Equivalent transfer function:');
disp(W_final);

[N,D] = numden(W_final);

Numerator = double(sym2poly(expand(N)));
Denominator = double(sym2poly(expand(D)));

fprintf('Numerator: ');
disp(Numerator);

fprintf('Denominator: ');
disp(Denominator);
