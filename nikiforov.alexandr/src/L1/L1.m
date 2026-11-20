clear; clc; close all;

K  = 0.7;
T  = 0.4;
T1 = 0.5;
T2 = 0.7;
T1_aper = 1.5; 

time = 0:0.005:5;
dt = time(2) - time(1);
s_tf = tf('s');

link_names = {
    '1. Интегрирующее'
    '2. Дифференцирующее'
    '3. Усилительное'
    '4. Апериодическое 1-го порядка'
    '5. Апериодическое 2-го порядка'
    '6. Колебательное'
    '7. Консервативное'
    '8. Инерционное интегрирующее'
    '9. Реальное дифференцирующее'
    '10. Форсирующее 1-го порядка'
    '11. Изодромное'
};

syms s positive;
syms t real;
syms K_sym T_sym T1_sym T2_sym positive;

W_sym = [
    K_sym / s
    K_sym * s
    K_sym
    K_sym / (T_sym * s + 1)
    K_sym / (T2_sym^2 * s^2 + T1_sym * s + 1)
    K_sym / (T2_sym^2 * s^2 + T1_sym * s + 1)
    K_sym / (T_sym^2 * s^2 + 1)
    K_sym / (s * (T_sym * s + 1))
    (K_sym * s) / (T_sym * s + 1)
    K_sym * (T_sym * s + 1)
    (K_sym * (T_sym * s + 1)) / s
];

U_step = laplace(heaviside(t), t, s);
Y_step = W_sym * U_step;
h = ilaplace(Y_step, s, t);

U_impulse = laplace(dirac(t), t, s);
Y_impulse = W_sym * U_impulse;
w = ilaplace(Y_impulse, s, t);

fprintf('Аналитические выражения h(t) и w(t)\n');
for i = 1:numel(W_sym)
    fprintf('%s:\n', link_names{i});
    fprintf('h(t) = %s\n', char(h(i)));
    fprintf('w(t) = %s\n\n', char(w(i)));
end

W_tf = {
    K / s_tf
    K * s_tf
    tf(K, 1)
    K / (T * s_tf + 1)
    K / (T2^2 * s_tf^2 + T1_aper * s_tf + 1)
    K / (T2^2 * s_tf^2 + T1 * s_tf + 1)
    K / (T^2 * s_tf^2 + 1)
    K / (s_tf * (T * s_tf + 1))
    (K * s_tf) / (T * s_tf + 1)
    K * (T * s_tf + 1)
    (K * (T * s_tf + 1)) / s_tf
};

figure('Name', 'Переходные характеристики h(t)');
for i = 1:numel(W_tf)
    subplot(4, 3, i);
    hold on; grid on;
    
    if i == 2
        y = zeros(size(time));
        y(1) = K / dt;
    elseif i == 3
        y = repmat(K, size(time));
    elseif i == 10
        y = repmat(K, size(time));
        y(1) = K + (K * T) / dt;
    else
        y = step(W_tf{i}, time);
    end
    
    plot(time, y, 'b', 'LineWidth', 1.2);
    title(link_names{i}, 'FontSize', 8);
    xlabel('t, с', 'FontSize', 7);
    ylabel('h(t)', 'FontSize', 7);
    xlim([0, 5]);
    hold off;
end

figure('Name', 'Импульсные характеристики w(t)');
for i = 1:numel(W_tf)
    subplot(4, 3, i);
    hold on; grid on;
    
    if i == 2
        y = zeros(size(time));
        y(1) = K / dt^2;
        y(2) = -K / dt^2;
    elseif i == 3
        y = zeros(size(time));
        y(1) = K / dt;
    elseif i == 10
        y = zeros(size(time));
        y(1) = (K * T) / dt^2;
        y(2) = -(K * T) / dt^2;
    else
        y = impulse(W_tf{i}, time);
    end
    
    plot(time, y, 'b', 'LineWidth', 1.2);
    title(link_names{i}, 'FontSize', 8);
    xlabel('t, с', 'FontSize', 7);
    ylabel('w(t)', 'FontSize', 7);
    xlim([0, 5]);
    hold off;
end

K_values = [0.3, 0.5, K, 1.0, 1.4];
time_long = 0:0.01:10;

figure('Name', 'Влияние коэффициента K');
hold on; grid on;
for k = K_values
    W_check = k / (T2^2 * s_tf^2 + T1 * s_tf + 1);
    step(W_check, time_long);
end
title('Колебательное звено: реакция при различных значениях K');
xlabel('t, с');
ylabel('h(t)');
legend('K = 0.3', 'K = 0.5', 'K = 0.7 (номинал)', 'K = 1.0', 'K = 1.4', 'Location', 'best');
hold off;

T_factors = [0.2, 0.5, 1.0, 2.0, 5.0];

figure('Name', 'Влияние постоянной времени T');
hold on; grid on;
for scale = T_factors
    W_check = K / (scale * T * s_tf + 1);
    step(W_check, time_long);
end
title('Апериодическое звено №4: масштабирование T в 2 и 5 раз');
xlabel('t, с');
ylabel('h(t)');
legend('T/5 (0.08 с)', 'T/2 (0.20 с)', 'T (0.40 с)', '2T (0.80 с)', '5T (2.00 с)', 'Location', 'best');
hold off;

T1_values = [1.2, 0.8, T1, 0.2, 0.05];
time_osc = 0:0.01:15;

figure('Name', 'Влияние коэффициента T1 на колебательное звено');
hold on; grid on;
for t1_val = T1_values
    W_check = K / (T2^2 * s_tf^2 + t1_val * s_tf + 1);
    step(W_check, time_osc);
end
title(sprintf('Колебательное звено №6: уменьшение T_1 при T_2 = %.1f с', T2));
xlabel('t, с');
ylabel('h(t)');
legend('T_1 = 1.2 с', 'T_1 = 0.8 с', 'T_1 = 0.5 с (номинал)', 'T_1 = 0.2 с', 'T_1 = 0.05 с', 'Location', 'best');
hold off;
