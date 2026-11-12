clear; clc; close all;

K_def = 1;
T = 1;
T1 = 0.2;
T2 = 0.7;
T1a = 2*T2;

time = 0:0.01:10;

names = {
    'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
    'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', 'Колебательное', ...
    'Консервативное', 'Инерционное интегрирующее', 'Реальное дифференцирующее', ...
    'Форсирующее', 'Изодромное'
};

N = numel(names);

syms s t K real
W_sym = {
    K/s, ...
    K*s, ...
    K, ...
    K/(T*s + 1), ...
    K/(T2^2*s^2 + T1a*s + 1), ...
    K/(T2^2*s^2 + T1*s + 1), ...
    K/(T*s^2 + 1), ...
    K/(s*(T*s + 1)), ...
    K*s/(T*s + 1), ...
    K*(T*s + 1), ...
    K*(T*s + 1)/s
};

U_step = laplace(heaviside(t),t,s);
U_impulse = laplace(dirac(t),t,s);

%% График переходных характеристик h(t)
figure('Name', 'Переходные характеристики h(t)', 'Color', 'w', 'Position', [50 50 1000 500]);
hold on; grid on;
for i = 1:N
    name = names{i};
    W_k = subs(W_sym{i}, K, K_def);
    Y_s = W_k * U_step;
    h_t = ilaplace(Y_s, s, t);
    y_vals = double(subs(h_t, t, time));
    y_vals(isinf(y_vals) | isnan(y_vals)) = NaN;
    plot(time, y_vals, 'DisplayName', name);
end
title('Переходные характеристики h(t)', 'FontSize', 14);
xlabel('t, с'); ylabel('h(t)'); xlim([0 10]);
legend('show', 'Location', 'bestoutside', 'FontSize', 9); hold off;

%% График импульсных характеристик w(t)
figure('Name', 'Импульсные характеристики w(t)', 'Color', 'w', 'Position', [50 50 1000 500]);
hold on; grid on;
for i = 1:N
    name = names{i};
    W_k = subs(W_sym{i}, K, K_def);
    Y_s = W_k * U_impulse;
    w_t = ilaplace(Y_s, s, t);
    y_vals = double(subs(w_t, t, time));
    y_vals(isinf(y_vals) | isnan(y_vals)) = NaN;
    plot(time, y_vals, 'DisplayName', name);
end
title('Импульсные характеристики w(t)', 'FontSize', 14);
xlabel('t, с'); ylabel('w(t)'); xlim([0 10]);
legend('show', 'Location', 'bestoutside', 'FontSize', 9); hold off;

%% Проверка с помощью Control System Toolbox
s_tf = tf('s');

W_num = {
    K_def/s_tf, ...
    K_def*s_tf, ...
    K_def, ...
    K_def/(T*s_tf + 1), ...
    K_def/(T2^2*s_tf^2 + T1a*s_tf + 1), ...
    K_def/(T2^2*s_tf^2 + T1*s_tf + 1), ...
    K_def/(T*s_tf^2 + 1), ...
    K_def/(s_tf*(T*s_tf + 1)), ...
    K_def*s_tf/(T*s_tf + 1), ...
    K_def*(T*s_tf + 1), ...
    K_def*(T*s_tf + 1)/s_tf
};

figure('Name', 'Проверка step()', 'Color', 'w', 'Position', [50 50 1000 500]);
hold on; grid on;
for i = 1:N
    try
        step(W_num{i}, time);
    catch ME
        warning('step() неприменимо для звена [%d] %s: %s', i, names{i}, ME.message);
    end
end
title('Проверка переходных характеристик', 'FontSize', 14);
xlabel('t, с'); ylabel('Выход'); xlim([0 10]); hold off;

figure('Name', 'Проверка impulse()', 'Color', 'w', 'Position', [50 50 1000 500]);
hold on; grid on;
for i = 1:N
    try
        impulse(W_num{i}, time);
    catch ME
        warning('impulse() неприменимо для звена [%d] %s: %s', i, names{i}, ME.message);
    end
end
title('Проверка импульсных характеристик', 'FontSize', 14);
xlabel('t, с'); ylabel('Выход'); xlim([0 10]); hold off;

%% * Влияние K на реальное дифференцирующее звено 
K_values = [0.5, 1.5, 2, 5];
idx_real_diff = 9;

figure('Name', 'Влияние K на реальное дифф. звено', 'Color', 'w', 'Position', [50 50 1000 500]);

subplot(1, 2, 1); hold on; grid on;
for k = K_values
    W_k = subs(W_sym{idx_real_diff}, K, k);
    y_h = double(subs(ilaplace(W_k * U_step, s, t), t, time));
    y_h(isinf(y_h)|isnan(y_h)) = NaN;
    plot(time, y_h, 'DisplayName', sprintf('K=%.1f', k));
end
title('Переходная функция h(t)', 'FontSize', 12);
xlabel('t, с'); ylabel('h(t)'); xlim([0 10]);
legend('show', 'Location', 'northeast', 'FontSize', 9);

subplot(1, 2, 2); hold on; grid on;
for k = K_values
    W_k = subs(W_sym{idx_real_diff}, K, k);
    y_w = double(subs(ilaplace(W_k * U_impulse, s, t), t, time));
    y_w(isinf(y_w)|isnan(y_w)) = NaN;
    plot(time, y_w, 'DisplayName', sprintf('K=%.1f', k));
end
title('Импульсная функция w(t)', 'FontSize', 12);
xlabel('t, с'); ylabel('w(t)'); xlim([0 10]);
legend('show', 'Location', 'northeast', 'FontSize', 9);

sgtitle('Реальное дифференцирующее звено: влияние K', 'FontSize', 14); hold off;

%% * Влияние увеличения и уменьшения T в 2 и в 5 раз на апериодическое звено 1-го порядка
T_values = [T/5, T/2, T, T*2, T*5];
labels_T = {'T/5', 'T/2', 'T', '2T', '5T'};

figure('Name', 'Влияние T на апериодическое звено', 'Color', 'w', 'Position', [50 50 1000 500]);
hold on; grid on;
for j = 1:numel(T_values)
    W_T = K_def/(T_values(j)*s_tf + 1);
    step(W_T, time);
    L = findobj(gca, 'Type', 'line');
    if ~isempty(L), L(1).DisplayName = labels_T{j}; end
end
title('Апериодическое звено 1-го порядка: влияние T', 'FontSize', 14);
xlabel('t, с'); ylabel('h(t)'); xlim([0 10]);
legend('show', 'Location', 'northeast', 'FontSize', 9); hold off;

%% * Влияние уменьшения T1 на колебательное звено
T1_values = [0.5, 0.1, 0.05, 0.01, 0.005];
T2_fix = T2;
labels_T1 = {'0.5', '0.1', '0.05', '0.01', '0.005'};

figure('Name', 'Влияние T1 на колебательное звено', 'Color', 'w', 'Position', [50 50 1000 500]);
hold on; grid on;
for j = 1:numel(T1_values)
    W_T1 = K_def/(T2_fix^2*s_tf^2 + T1_values(j)*s_tf + 1);
    step(W_T1, time);
    L = findobj(gca, 'Type', 'line');
    if ~isempty(L), L(1).DisplayName = labels_T1{j}; end
end
title('Колебательное звено: влияние T1 при фиксированном T2');
xlabel('t, с'); ylabel('h(t)'); xlim([0 10]);
legend('show', 'Location', 'northeast', 'FontSize', 9); hold off;
