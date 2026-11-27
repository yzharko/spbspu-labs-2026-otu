clear;
clc;
close all;
format short g;

%% 1. Исходные данные
n = 2;
A = [0, 2; -10, -2];
B = [1; 1];
CT = [5; 0];
C = CT.';
D = 0;
M = [0.5, 0; 6, 2];
tolerance = 1e-9;

%% 2. Замена базиса x = M*x_new
A_new = M \ (A*M);
B_new = M \ B;
C_new = C*M;

fprintf('det(M) = %g\n', det(M));
disp('A_new ='); disp(A_new);
disp('B_new ='); disp(B_new);
disp('C_new ='); disp(C_new);

%% 3. Передаточная функция
den_analytic = [1, -trace(A), det(A)];
num_analytic = [0, C*B, C*(A - trace(A)*eye(n))*B];
W_tf = tf(num_analytic, den_analytic);

[num, den] = ss2tf(A, B, C, D);
[num_new, den_new] = ss2tf(A_new, B_new, C_new, D);
assert(norm(num-num_analytic, inf) < tolerance && ...
       norm(den-den_analytic, inf) < tolerance, ...
       'Аналитические коэффициенты расходятся с ss2tf.');
assert(norm(num_new-num, inf) < tolerance && ...
       norm(den_new-den, inf) < tolerance, ...
       'Передаточная функция после замены базиса изменилась.');

fprintf('W(s) = (%g*s + %g)/(s^2 + %g*s + %g)\n', ...
    num_analytic(2), num_analytic(3), den_analytic(2), den_analytic(3));
fprintf('ss2tf, исходный базис: num = %s, den = %s\n', ...
    mat2str(num, 12), mat2str(den, 12));
fprintf('ss2tf, новый базис:    num = %s, den = %s\n', ...
    mat2str(num_new, 12), mat2str(den_new, 12));
disp('Полюса W(s):'); disp(pole(W_tf));
fprintf('W(0) = %.12g\n', dcgain(W_tf));

%% 4. Переходные характеристики тремя способами
t = (0:0.01:10).';
sys_old = ss(A, B, C, D);
sys_new = ss(A_new, B_new, C_new, D);
y_old = step(sys_old, t);
y_new = step(sys_new, t);
y_tf = step(W_tf, t);

fig_main = figure('Name', 'Сравнение трёх моделей', 'Color', 'w', ...
                  'Position', [120, 120, 840, 410]);
plot(t, y_old, '-', 'Color', [0.08, 0.25, 0.49], 'LineWidth', 2.6);
hold on;
plot(t, y_new, '--', 'Color', [0.78, 0.26, 0.09], 'LineWidth', 1.6);
plot(t, y_tf, 'o', 'Color', [0.10, 0.48, 0.31], 'LineWidth', 1, ...
    'MarkerIndices', 1:24:numel(t), 'MarkerSize', 4.3);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k');
grid on;
xlabel('t, с'); ylabel('y(t)');
xlim([0, 10]); ylim([0, 2]);
legend('Исходная модель', 'Модель в новом базисе', 'Передаточная функция', ...
       'Location', 'southeast', 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.4, 0.4, 0.4]);

err_new = max(abs(y_old-y_new));
err_tf = max(abs(y_old-y_tf));
fprintf('max |y - y_new|  = %.6e\n', err_new);
fprintf('max |y - y_tf|   = %.6e\n', err_tf);
assert(max([err_new, err_tf]) < tolerance, ...
    'Переходные характеристики трёх моделей расходятся.');

%% 5. Критерий Калмана
U = ctrb(A, B);
N = obsv(A, C);
rank_U = rank(U);
rank_N = rank(N);
disp('Матрица управляемости U ='); disp(U);
fprintf('det(U) = %g, rank(U) = %d\n', det(U), rank_U);
disp('Матрица наблюдаемости N ='); disp(N);
fprintf('det(N) = %g, rank(N) = %d\n', det(N), rank_N);
assert(rank_U == n && rank_N == n, ...
    'Для заданных преобразований нужны полные ранги U и N.');
disp('Система полностью управляема и полностью наблюдаема.');

%% 6. Канонические формы
a1 = den_analytic(2);
T_c = [A*B + a1*B, B];
A_c = T_c \ (A*T_c);
B_c = T_c \ B;
C_c = C*T_c;

% Наблюдаемая форма двойственна управляемой
A_o_expected = A_c.';
B_o_expected = C_c.';
C_o_expected = B_c.';
T_o = N \ obsv(A_o_expected, C_o_expected);
A_o = T_o \ (A*T_o);
B_o = T_o \ B;
C_o = C*T_o;

A_c_expected = [0, 1; -den_analytic(3), -a1];
assert(norm(A_c - A_c_expected, inf) < tolerance && ...
       norm(B_c - [0; 1], inf) < tolerance && ...
       norm(C_c - fliplr(num_analytic(2:end)), inf) < tolerance && ...
       norm(A_o - A_o_expected, inf) < tolerance && ...
       norm(B_o - B_o_expected, inf) < tolerance && ...
       norm(C_o - C_o_expected, inf) < tolerance, ...
       'Матрицы преобразованных моделей не совпадают с каноническими.');

disp('T_c ='); disp(T_c);
disp('A_c ='); disp(A_c);
disp('B_c ='); disp(B_c);
disp('C_c ='); disp(C_c);
disp('T_o ='); disp(T_o);
disp('A_o ='); disp(A_o);
disp('B_o ='); disp(B_o);
disp('C_o ='); disp(C_o);

%% 7. Сравнение откликов канонических моделей
y_c = step(ss(A_c, B_c, C_c, D), t);
y_o = step(ss(A_o, B_o, C_o, D), t);

fig_canonical = figure('Name', 'Сравнение канонических форм', 'Color', 'w', ...
                       'Position', [160, 150, 840, 410]);
plot(t, y_old, '-', 'Color', [0.08, 0.25, 0.49], 'LineWidth', 2.6);
hold on;
plot(t, y_c, '--', 'Color', [0.78, 0.26, 0.09], 'LineWidth', 1.6);
plot(t, y_o, 'o', 'Color', [0.10, 0.48, 0.31], 'LineWidth', 1, ...
    'MarkerIndices', 1:24:numel(t), 'MarkerSize', 4.3);
set(gca, 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'GridColor', 'k');
grid on;
xlabel('t, с'); ylabel('y(t)');
xlim([0, 10]); ylim([0, 2]);
legend('Исходная модель', 'Управляемая форма', 'Наблюдаемая форма', ...
       'Location', 'southeast', 'Color', 'w', 'TextColor', 'k', 'EdgeColor', [0.4, 0.4, 0.4]);

err_c = max(abs(y_old-y_c));
err_o = max(abs(y_old-y_o));
fprintf('max |y - y_c|    = %.6e\n', err_c);
fprintf('max |y - y_o|    = %.6e\n', err_o);
assert(max([err_c, err_o]) < tolerance, ...
    'Переходные характеристики канонических моделей расходятся.');
