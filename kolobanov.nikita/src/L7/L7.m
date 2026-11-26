clear;
clc;
close all;

A = [-0.5 2; -12 -4];
B = [1; 1];
C = [5 0.5];
D = 0;
M = [1 5; 0 2];
ORDER = size(A, 1);
TOLERANCE = 1e-9;
TIME = (0:0.005:4).';
MARKER_INDICES = 1:20:numel(TIME);

A_new = M \ (A*M);
B_new = M \ B;
C_new = C*M;
disp('Новый базис [A B; C D]');
disp([A_new B_new; C_new D]);

syms s;
E = eye(ORDER);
sE_A = s*E - A;
W_sym = simplify(C * (sE_A \ B));
fprintf('W(s) аналитически = %s\n', char(W_sym));
[num, den] = ss2tf(A, B, C, D);
[num_new, den_new] = ss2tf(A_new, B_new, C_new, D);
W_tf = tf(num, den);
[sym_num, sym_den] = numden(W_sym);
sym_den = sym2poly(sym_den);
analytic_coeff_error = max(abs([num den] - [0 sym2poly(sym_num) sym_den]/sym_den(1)));
basis_coeff_error = max(abs([num den] - [num_new den_new]));
disp(table([num; den], [num_new; den_new], ...
    'VariableNames', {'OriginalBasis', 'NewBasis'}, ...
    'RowNames', {'Numerator', 'Denominator'}));
fprintf('Расхождения: аналитика %.4g, базисы %.4g. Допуск %.0e соблюдён: %s.\n', ...
    analytic_coeff_error, basis_coeff_error, TOLERANCE, ...
    mat2str(max(analytic_coeff_error, basis_coeff_error) < TOLERANCE));

sys_old = ss(A, B, C, D);
sys_new = ss(A_new, B_new, C_new, D);
y_old = step(sys_old, TIME);
y_new = step(sys_new, TIME);
y_tf = step(W_tf, TIME);
figure('Name', 'Сравнение трёх моделей', 'Color', 'w');
plot(TIME, y_old, 'k-', 'LineWidth', 3);
hold on;
plot(TIME, y_new, 'b--o', 'LineWidth', 2, 'MarkerIndices', MARKER_INDICES(1:2:end));
plot(TIME, y_tf, 'r:^', 'LineWidth', 1.25, 'MarkerIndices', MARKER_INDICES(2:2:end));
grid on;
xlabel('t, с');
ylabel('Выход y(t), отн. ед.');
legend('Исходная модель', 'Новый базис', 'Передаточная функция', ...
    'Location', 'southoutside', 'Orientation', 'horizontal');

U = ctrb(A, B);
N = obsv(A, C);
disp(table(U, N));
fprintf('Ранг матрицы управляемости %d из %d.\n', rank(U), ORDER);
fprintf('Ранг матрицы наблюдаемости %d из %d.\n', rank(N), ORDER);

a1 = den(2);
b1 = num(2);
Ac = [0 -den(3); 1 -a1];
Bc = [1; 0];
Cc = [b1 num(3)-a1*b1];
controllable = ss(Ac, Bc, Cc, D);
observable = ss(Ac.', Cc.', Bc.', D);
disp('[A B; C D] (управляемая, наблюдаемая)');
disp([Ac Bc; Cc D; observable.A observable.B; observable.C D]);

y_control = step(controllable, TIME);
y_observe = step(observable, TIME);
max_error = max(abs([y_new y_tf y_control y_observe] - y_old)).';
disp(table(max_error, 'RowNames', {'NewBasis', 'Transfer', 'Controllable', 'Observable'}));
figure('Name', 'Сравнение канонических форм', 'Color', 'w');
plot(TIME, y_old, 'k-', 'LineWidth', 3);
hold on;
plot(TIME, y_control, 'b--o', 'LineWidth', 2, 'MarkerIndices', MARKER_INDICES(1:2:end));
plot(TIME, y_observe, 'r:^', 'LineWidth', 1.25, 'MarkerIndices', MARKER_INDICES(2:2:end));
grid on;
xlabel('t, с');
ylabel('Выход y(t), отн. ед.');
legend('Исходная модель', 'Управляемая форма', 'Наблюдаемая форма', ...
    'Location', 'southoutside', 'Orientation', 'horizontal');
