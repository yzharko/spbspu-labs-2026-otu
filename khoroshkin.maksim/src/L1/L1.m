clear; clc; close all;

syms s K T positive
syms t real

K_var = 0.9;
T_val = 0.6;

% T1 >= 2*T2
T1_5 = 2;   T2_5 = 0.7;

% T1 < 2*T2
T1_6 = 0.8; T2_6 = 0.7;

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

W = {K/s, K*s, K, ...
    K/(T*s + 1), ...
    K/(T2_5^2*s^2 + T1_5*s + 1), ...
    K/(T2_6^2*s^2 + T1_6*s + 1), ...
    K/(T*s^2 + 1), ...
    K/(s*(T*s + 1)), ...
    K*s/(T*s + 1), ...
    K*(T*s + 1), ...
    K*(T*s + 1)/s};

h = cell(1, numel(W));
w = cell(1, numel(W));
for i = 1:numel(W)
    h{i} = ilaplace(W{i}*U_step, s, t);
    w{i} = ilaplace(W{i}*U_impulse, s, t);
end

disp('h(t):');
disp(h.');
disp('w(t):');
disp(w.');

h_plot = cell(1, numel(W));
w_plot = cell(1, numel(W));
for i = 1:numel(W)
    h_plot{i} = subs(h{i}, {dirac(t), dirac(1, t)}, {0, 0});
    w_plot{i} = subs(w{i}, {dirac(t), dirac(1, t)}, {0, 0});
end

names = {'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
    'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', ...
    'Колебательное', 'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее 1-го порядка', 'Изодромное'};
n = numel(names);
t_end = 10;
time = 0:0.01:t_end;

h_label = {'', 'h(t) = K\delta(t)', '', '', '', '', '', '', '', ...
    'в нуле: KT\delta(t)', ''};
w_label = {'', 'w(t) = K\delta''(t)', 'w(t) = K\delta(t)', '', '', '', '', '', ...
    'в нуле: (K/T)\delta(t)', 'w(t) = K\delta(t) + KT\delta''(t)', ...
    'в нуле: KT\delta(t)'};

figure('Name', 'Переходные функции h(t)');
for i = 1:n
    subplot(3, 4, i);
    hold on;
    grid on;
    h_k = subs(h_plot{i}, [K T], [K_var T_val]);
    y = double(subs(h_k, t, time));
    plot(time, y);
    if ~isempty(h_label{i})
        text(0.05, 0.85, h_label{i}, 'Units', 'normalized', 'FontSize', 9);
    end
    title(sprintf('%d. %s', i, names{i}), 'FontSize', 9);
    xlabel('Время, с');
    ylabel('h(t)');
    hold off;
end

figure('Name', 'Функции веса w(t)');
for i = 1:n
    subplot(3, 4, i);
    hold on;
    grid on;
    w_k = subs(w_plot{i}, [K T], [K_var T_val]);
    y = double(subs(h_k, t, time));
    plot(time, y);
    if ~isempty(w_label{i})
        text(0.05, 0.85, w_label{i}, 'Units', 'normalized', 'FontSize', 9);
    end
    title(sprintf('%d. %s', i, names{i}), 'FontSize', 9);
    xlabel('Время, с');
    ylabel('w(t)');
    hold off;
end

s_tf = tf('s');
W_tf = {K_var/s_tf, [], tf(K_var, 1), ...
    K_var/(T_val*s_tf + 1), ...
    K_var/(T2_5^2*s_tf^2 + T1_5*s_tf + 1), ...
    K_var/(T2_6^2*s_tf^2 + T1_6*s_tf + 1), ...
    K_var/(T_val*s_tf^2 + 1), ...
    K_var/(s_tf*(T_val*s_tf + 1)), ...
    K_var*s_tf/(T_val*s_tf + 1), ...
    [], ...
    K_var*(T_val*s_tf + 1)/s_tf};

figure('Name', 'Проверка: step');
for i = 1:n
    if ~isempty(W_tf{i})
        subplot(3, 4, i);
        y = step(W_tf{i}, time);
        plot(time, y);
        if max(y) - min(y) < 1e-6
            ylim([min(y) - 1, max(y) + 1]);
        end
        grid on;
        title(names{i});
    end
end

figure('Name', 'Проверка: impulse');
for i = 1:n
    if ~isempty(W_tf{i}) && i ~= 3
        subplot(3, 4, i);
        y = impulse(W_tf{i}, time);
        plot(time, y);
        if max(y) - min(y) < 1e-6
            ylim([min(y) - 1, max(y) + 1]);
        end
        grid on;
        title(names{i});
    end
end

T_factors = [1/5 1/2 1 2 5];

figure('Name', 'Апериодическое звено: вариация T');
hold on;
grid on;
fprintf('Апериодическое звено, K = %.1f:\n', K_var);
for m = T_factors
    T_m = T_val * m;
    sys = tf(K_var, [T_m 1]);
    info = stepinfo(sys, 'SettlingTimeThreshold', 0.05);
    fprintf('T = %5.2f   время установления = %6.2f   3T = %6.2f\n', ...
        T_m, info.SettlingTime, 3*T_m);
    [y, tt] = step(sys, time);
    plot(tt, y);
end
legend(compose('T = %.2f', T_val * T_factors), 'Location', 'southeast');
title('Апериодическое звено: влияние T');
xlabel('Время, с');
ylabel('h(t)');
hold off;

T1_list = [1.2 0.8 0.4 0.2 0.1];

figure('Name', 'Колебательное звено: вариация T1');
hold on;
grid on;
fprintf('\nКолебательное звено, K = %.1f, T2 = %.1f:\n', K_var, T2_6);
for T1_cur = T1_list
    sys = tf(K_var, [T2_6^2 T1_cur 1]);
    info = stepinfo(sys, 'SettlingTimeThreshold', 0.05);
    zeta = T1_cur / (2*T2_6);
    sigma_theory = 100 * exp(-pi*zeta / sqrt(1 - zeta^2));
    fprintf('T1 = %4.2f  zeta = %.3f  перерегулирование = %5.1f%%  время затухания = %6.1f\n', ...
        T1_cur, zeta, info.Overshoot, sigma_theory, info.SettlingTime);
    [y, tt] = step(sys, 0:0.01:40);
    plot(tt, y);
end
legend(compose('T1 = %.1f', T1_list), 'Location', 'northeast');
title('Колебательное звено: влияние T_1 при T_2 = 0.7');
xlabel('Время, с');
ylabel('h(t)');
hold off;

K_list = [0.5 1 2 4];

figure('Name', 'Колебательное звено: вариация K');
hold on;
grid on;
for K_cur = K_list
    [y, tt] = step(tf(K_cur, [T2_6^2 T1_6 1]), time);
    plot(tt, y);
end
legend(compose('K = %.1f', K_list), 'Location', 'northeast');
title('Колебательное звено: влияние K');
xlabel('Время, с');
ylabel('h(t)');
hold off;
