clear;
clc;
close all;
A  = [-2,   4;
       0.5, -4];
B  = [1; 1];
CT = [1; 0];
C  = CT';
M  = [2, 1;
      0, 4];
n = size(A, 1);
det_M = det(M);
if abs(det_M) < 1e-12
    error('Матрица преобразования M вырождена (det(M) = 0).');
end

M_inv = inv(M);
A_new = M_inv * A * M;
B_new = M_inv * B;
C_new = C * M;
fprintf('det(M) = %.6f\n', det_M);
fprintf('A_new =\n'); disp(A_new);
fprintf('B_new =\n'); disp(B_new);
fprintf('C_new =\n'); disp(C_new);

syms s;
sI_A = s * eye(n) - A;
W_sym = simplify(C * inv(sI_A) * B);
fprintf('W(s) аналитически = %s\n', char(W_sym));

[num_old, den_old] = ss2tf(A, B, C, 0); %численно переводит исходную модель пространства состояний в коэффициенты числителя и знаменателя передаточной функции
W_tf_old = tf(num_old, den_old);
disp('W(s) по исходной модели (tf):'); disp(W_tf_old);

[num_new, den_new] = ss2tf(A_new, B_new, C_new, 0); 
W_tf_new = tf(num_new, den_new);
disp('W_new(s) по модели в новом базисе (tf):'); disp(W_tf_new);

sys_old = ss(A, B, C, 0);
sys_new = ss(A_new, B_new, C_new, 0);
p = pole(sys_old);
max_tau = max(1 ./ abs(real(p(real(p) < 0)))); %Считает постоянные времени и выбирает наибольшую, признак затухания системы
t_end = max(5 * max_tau, 1);
t = 0:(t_end / 1000):t_end;

[y_old, ~] = step(sys_old, t);
[y_new, ~] = step(sys_new, t);
[y_tf,  ~] = step(W_tf_old, t);

figure('Name', 'Переходные процессы', 'Color', 'w');
plot(t, y_old, 'r-',  'LineWidth', 1.5); hold on;
plot(t, y_new, 'b--', 'LineWidth', 1.5);
plot(t, y_tf,  'g:',  'LineWidth', 2.0);
grid on;
xlabel('Время t, с');
ylabel('Выходной сигнал y(t)');
title('Переходный процесс при единичном ступенчатом воздействии');
legend('Исходная модель (A, B, C)', ...
       'Модель в новом базисе (A_{нов}, B_{нов}, C_{нов})', ...
       'Передаточная функция W(s)', ...
       'Location', 'southeast');

U = ctrb(A, B);
rank_U = rank(U);
fprintf('\nМатрица управляемости U =\n'); disp(U);
fprintf('det(U) = %.6f, rank(U) = %d при n = %d\n', det(U), rank_U, n);
N = obsv(A, C);
rank_N = rank(N);
fprintf('Матрица наблюдаемости N =\n'); disp(N);
fprintf('det(N) = %.6f, rank(N) = %d при n = %d\n', det(N), rank_N, n);

a1 = den_old(2);
a0 = den_old(3);
b1 = num_old(2);
b0 = num_old(3);
A_c_target = [0, 1; -a0, -a1];
B_c_target = [0; 1];
U_c = ctrb(A_c_target, B_c_target);
P_c = U / U_c;               
P_c_inv = inv(P_c);
A_c = P_c_inv * A * P_c;
B_c = P_c_inv * B;
C_c = C * P_c;
fprintf('\n Каноническая управляемая форма \n');
fprintf('Матрица преобразования P_c =\n'); disp(P_c);
fprintf('A_c =\n'); disp(A_c);
fprintf('B_c =\n'); disp(B_c);
fprintf('C_c =\n'); disp(C_c);

A_o_target = [0, -a0; 1, -a1];
C_o_target = [0, 1];
N_o = obsv(A_o_target, C_o_target);
P_o = N \ N_o;               
P_o_inv = inv(P_o);
A_o = P_o_inv * A * P_o;
B_o = P_o_inv * B;
C_o = C * P_o;
fprintf('\n Каноническая наблюдаемая форма \n');
fprintf('Матрица преобразования P_o =\n'); disp(P_o);
fprintf('A_o =\n'); disp(A_o);
fprintf('B_o =\n'); disp(B_o);
fprintf('C_o =\n'); disp(C_o);

%Моделирование и сопоставление переходных процессов канонических форм
sys_c = ss(A_c, B_c, C_c, 0);
sys_o = ss(A_o, B_o, C_o, 0);
[y_c, ~] = step(sys_c, t);
[y_o, ~] = step(sys_o, t);
max_diff_canon = max([max(abs(y_old - y_c)), max(abs(y_old - y_o))]);
fprintf('Максимальное расхождение переходных процессов канонических форм: %.2e\n', max_diff_canon);
figure('Name', 'Сопоставление канонических форм', 'Color', 'w');
plot(t, y_old, 'r-',  'LineWidth', 1.5); hold on;
plot(t, y_c,   'b--', 'LineWidth', 1.5);
plot(t, y_o,   'k:',  'LineWidth', 2.0);
grid on;
xlabel('Время t, с');
ylabel('Выходной сигнал y(t)');
title('Сравнение переходных характеристик исходной модели, КУФ и КНФ');
legend('Исходная модель', 'Каноническая управляемая форма', 'Каноническая наблюдаемая форма', ...
       'Location', 'southeast');
