clear; clc; close all;

% параметры варианта
K  = 3.2;
T  = 0.9;
T1 = 1.0;
T2 = 1.4;
T1_ap = 2 * T2; % для апериодического 2-го порядка
eps_val = 0.001;
time = 0:0.01:10;
t_end = time(end);

syms s t real
assumeAlso(t >= 0);

names = {
    'Интегрирующее';
    'Идеальное дифференцирующее';
    'Усилительное';
    'Апериодическое 1го порядка';
    'Апериодическое 2го порядка';
    'Колебательное';
    'Консервативное';
    'Инерционное интегрирующее';
    'Реальное дифференцирующее';
    'Форсирующее 1го порядка';
    'Изодромное'};

W_sym = {
    K/s;
    K*s;
    K;
    K/(T*s + 1);
    K/(T2^2*s^2 + T1_ap*s + 1);
    K/(T2^2*s^2 + T1*s + 1);
    K/(T*s^2 + 1);
    K/(s*(T*s + 1));
    K*s/(T*s + 1);
    K*(T*s + 1);
    K*(T*s + 1)/s};

W = cell(1, 11);
W{1}  = tf(K, [1 0]);
W{2}  = tf([K 0], [eps_val 1]);
W{3}  = tf(K, 1);
W{4}  = tf(K, [T 1]);
W{5}  = tf(K, [T2^2, T1_ap, 1]);
W{6}  = tf(K, [T2^2, T1, 1]);
W{7}  = tf(K, [T 0 1]);
W{8}  = tf(K, [T 1 0]);
W{9}  = tf([K 0], [T 1]);
W{10} = tf([K*T K], [eps_val 1]);
W{11} = tf([K*T K], [1 0]);

U_step    = laplace(sym(1),   t, s);
U_impulse = laplace(dirac(t), t, s);

h = cell(1, numel(W_sym));
w = cell(1, numel(W_sym));
for n = 1:numel(W_sym)
    h{n} = simplify(ilaplace(W_sym{n} * U_step,    s, t));
    w{n} = simplify(ilaplace(W_sym{n} * U_impulse, s, t));
end

plot_grid(h, names, 'Переходные характеристики h(t)', t_end, 'h(t)');
plot_grid(w, names, 'Импульсные характеристики w(t)', t_end, 'w(t)');

% проверка
fprintf('\nПроверка результатов!!!\n');
for n = 1:numel(W_sym)
    h_func = matlabFunction(h{n}, 'Vars', t);
    y_anal = h_func(time);
    if isscalar(y_anal), y_anal = repmat(y_anal, size(time)); end
    y_cst = step(W{n}, time);
    err = max(abs(y_anal(:) - y_cst(:)));
    fprintf('Звено №%-2d (%-30s): погрешность = %.2e\n', n, names{n}, err);
end

% влияние T на звено №4
figure('Name', 'Влияние T на звено №4');
hold on; grid on;
T_factors = [1/5, 1/2, 1, 2, 5];
for f = T_factors
    step(tf(K, [T*f 1]), time);
end
title('Апериодическое звено №4: изменение T в 2 и 5 раз');
xlabel('Время t, с'); ylabel('h(t)');
legend('T/5', 'T/2', 'T', '2T', '5T');
hold off;

% влияние T1 на звено №6
figure('Name', 'Влияние T1 на звено №6');
hold on; grid on;
T1_values = [T1, 0.5, 0.2, 0.05];
for t1v = T1_values
    step(tf(K, [T2^2, t1v, 1]), time);
end
title('Колебательное звено №6: приближение T1 к 0');
xlabel('Время t, с'); ylabel('h(t)');
legend(compose('T1 = %.2f', T1_values));
hold off;

% замена всех сингулярных слагаемых на ноль
function expr = strip_singular(expr, t)
    expr = subs(expr, dirac(t), 0);
    for k = 1:2
        expr = subs(expr, diff(dirac(t), t, k), 0);
    end
end

function plot_grid(cells, names, fig_name, t_end, ylab)
    figure('Name', fig_name, 'Color', 'w');
    tiledlayout(4, 3, 'TileSpacing', 'compact');
    tt = linspace(0, t_end, 1000);
    for n = 1:numel(cells)
        nexttile;
        expr = strip_singular(cells{n}, sym('t'));
        y = double(subs(expr, sym('t'), tt));
        plot(tt, y, 'LineWidth', 1.4);
        grid on;
        title(sprintf('%d. %s', n, names{n}), 'Interpreter', 'none');
        xlabel('t, с'); ylabel(ylab);
    end
    sgtitle(fig_name);
end
