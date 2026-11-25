clear;
clc;
close all;

%% 1. Исходные данные
n = 2;
A  = [7, 1;
      3, 2];
B  = [6; 7];
CT = [6; 9];
C  = CT';      
M  = [8, 7;
      1, 4];

%% 2. Замена базиса
M_inv = inv(M);
A_new = M_inv * A * M;
B_new = M_inv * B;
C_new = C * M;

%% 3. Переход к передаточной функции
syms s;
sI_A = s*eye(n) - A;
W_sym = simplify(C * inv(sI_A) * B);
fprintf('W(s) аналитически = %s\n', char(W_sym));

[num, den] = ss2tf(A, B, C, 0);
W_tf = tf(num, den);
disp('W(s) (tf) = '); disp(W_tf);

[ya_num, ya_dem] = ss2tf(A_new, B_new, C_new, 0);
ya_W_tf = tf(ya_num, ya_dem);
disp('ya_W(s) (tf) = '); disp(ya_W_tf);
fprintf('Передаточная функция новой модели совпадает: %s\n', ...
    mat2str(all(abs(num-ya_num) < 1e-9) && all(abs(den-ya_dem) < 1e-9)));

%% 4. Моделирование переходного процесса
t = 0:0.01:10;
sys_old = ss(A, B, C, 0);
sys_new = ss(A_new, B_new, C_new, 0);

[y_old, ~] = step(sys_old, t);
[y_new, ~] = step(sys_new, t);
[y_tf, ~] = step(W_tf, t);

figure('Name','Переходные процессы (3 способа)','Color','w');
plot(t, y_old, 'r-', 'LineWidth', 1.5); hold on;
plot(t, y_new, 'b--', 'LineWidth', 1.5);
plot(t, y_tf, 'g:', 'LineWidth', 2);
grid on; xlabel('t, с'); ylabel('y(t)');
title('Переходной процесс при единичном ступенчатом воздействии');
legend('По исходной модели (A,B,C)','По модели в новом базисе','По передаточной функции W(s)','Location','southeast');

%% 5. Управляемость и Наблюдаемость
U = ctrb(A, B);
rank_U = rank(U);
fprintf('\nМатрица управляемости Uc =\n'); disp(U);
fprintf('rank(U) = %d при n = %d -> система %s\n', rank_U, n, ...
    char(ternary(rank_U == n, "управляема", "НЕ управляема")));

N = obsv(A, C);
rank_N = rank(N);
fprintf('\nМатрица наблюдаемости Uo =\n'); disp(N);
fprintf('rank(N) = %d при n = %d -> система %s\n', rank_N, n, ...
    char(ternary(rank_N == n, "наблюдаема", "НЕ наблюдаема")));

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
