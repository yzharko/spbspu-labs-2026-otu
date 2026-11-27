clear; close all; clc;

s   = tf('s');
W1  = 1 / (s + 1);             % W1(s) = 1/(s+1)
W2  = tf(2);                   % W2(s) = 2
W3  = 5 / (s^2 + 3*s + 2);     % W3(s) = 5/(s^2+3s+2)

t_end = 20;                    % время моделирования, с
tv    = linspace(0, t_end, 4000);

syms s_sym U X Y t
W1s = 1/(s_sym + 1);
W2s = sym(2);
W3s = 5/(s_sym^2 + 3*s_sym + 2);

sol = solve([X == U - Y, Y == W3s*(W2s*(W1s*X - Y) + X - Y)], [X, Y]);
Weq_sym = simplify(sol.Y / U);

[Ns, Ds] = numden(expand(Weq_sym));
fprintf('Символьное решение уравнений схемы:\n');
fprintf('  числитель   W_eq(s): %s\n', char(Ns));
fprintf('  знаменатель W_eq(s): %s\n\n', char(Ds));

W3i = feedback(W3, 1);


F = series(W3i, parallel(series(W1, W2), 1));   % W3i*(1 + W1*W2)
G = feedback(F, feedback(W2, W1));              % G = F/(1 + F*W2/(1+W1*W2))

Weq = minreal(feedback(G, 1));

Weq_hand = tf([5 15], [1 4 25 32]);

fprintf('ПФ, полученная функциями series/parallel/feedback:\n');
Weq
fprintf('ПФ, полученная вручную:\n');
Weq_hand

% Сравнение на сетке частот
sv = [0, 0.5, 1, 2, 5, 10];
err = 0;
for w = sv
    err = max(err, abs(evalfr(Weq, 1i*w) - evalfr(Weq_hand, 1i*w)));
end
fprintf('Максимальное отклонение ПФ на частотной сетке: %.2e\n\n', err);

fprintf('Полюса эквивалентной ПФ:\n');
disp(pole(Weq).');

Uin    = {1/s_sym, sym(1), 2/(s_sym^2 + 4)};    % изображения воздействий
lbl_in = {'Единичный скачок', 'Единичный импульс', 'Синусоида sin(2t)'};

figure('Color', 'w');
for k = 1:3
    Yk    = simplify(subs(sol.Y, U, Uin{k}));    % реакция исходной схемы
    y_or  = matlabFunction(ilaplace(Yk, s_sym, t), 'Vars', t);
    switch k
        case 1, y_eq = step(Weq, tv);
        case 2, y_eq = impulse(Weq, tv);
        otherwise, y_eq = lsim(Weq, sin(2*tv), tv);
    end
    subplot(3, 1, k); hold on; grid on;
    plot(tv, y_eq, 'r--', 'LineWidth', 1.5);
    yo = y_or(tv);
    plot(tv, yo, 'b-',  'LineWidth', 1.0);
    title(lbl_in{k}); xlabel('t, с'); ylabel('y(t)');
    legend('Эквивалентная ПФ (step/impulse/lsim)', ...
           'Исходная схема (символьное решение)', 'Location', 'best');
    fprintf('%s: max|y_исх - y_экв| = %.2e\n', ...
            lbl_in{k}, max(abs(y_eq(:) - yo(:))));
end
