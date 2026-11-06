clear; clc; close all;

K_val  = 2.6;
T_val  = 0.3;
T1_val = 0.6;
T2_val = 0.5;

T1_aper_val = 1.2; 

time = 0:0.005:10;
s_tf = tf('s');

syms s t positive;
syms K_sym T_sym T1_sym T2_sym positive;

U_step_sym = laplace(heaviside(t), t, s);

U_impulse_sym = laplace(dirac(t), t, s);

links = {
    struct('id', 1, 'name', 'Интегрирующее (идеальное)', ...
           'sym_W', K_sym / s, ...
           'tf_val', K_val / s_tf, ...
           'has_ideal_impulse', false), ...
    
    struct('id', 2, 'name', 'Дифференцирующее (идеальное)', ...
           'sym_W', K_sym * s, ...
           'tf_val', K_val * s_tf / (1e-4 * s_tf + 1), ...
           'has_ideal_impulse', true), ...
           
    struct('id', 3, 'name', 'Усилительное (безынерционное)', ...
           'sym_W', K_sym, ...
           'tf_val', tf(K_val, 1), ...
           'has_ideal_impulse', true), ...
           
    struct('id', 4, 'name', 'Апериодическое 1-го порядка', ...
           'sym_W', K_sym / (T_sym * s + 1), ...
           'tf_val', K_val / (T_val * s_tf + 1), ...
           'has_ideal_impulse', false), ...
           
    struct('id', 5, 'name', 'Апериодическое 2-го порядка (T1 >= 2*T2)', ...
           'sym_W', K_sym / (T2_sym^2 * s^2 + T1_sym * s + 1), ...
           'tf_val', K_val / (T2_val^2 * s_tf^2 + T1_aper_val * s_tf + 1), ...
           'has_ideal_impulse', false), ...
           
    struct('id', 6, 'name', 'Колебательное (T1 < 2*T2)', ...
           'sym_W', K_sym / (T2_sym^2 * s^2 + T1_sym * s + 1), ...
           'tf_val', K_val / (T2_val^2 * s_tf^2 + T1_val * s_tf + 1), ...
           'has_ideal_impulse', false), ...
           
    struct('id', 7, 'name', 'Консервативное', ...
           'sym_W', K_sym / (T_sym^2 * s^2 + 1), ...
           'tf_val', K_val / (T_val^2 * s_tf^2 + 1), ...
           'has_ideal_impulse', false), ...
           
    struct('id', 8, 'name', 'Инерционное интегрирующее', ...
           'sym_W', K_sym / (s * (T_sym * s + 1)), ...
           'tf_val', K_val / (s_tf * (T_val * s_tf + 1)), ...
           'has_ideal_impulse', false), ...
           
    struct('id', 9, 'name', 'Реальное дифференцирующее', ...
           'sym_W', (K_sym * s) / (T_sym * s + 1), ...
           'tf_val', (K_val * s_tf) / (T_val * s_tf + 1), ...
           'has_ideal_impulse', false), ...
           
    struct('id', 10, 'name', 'Форсирующее 1-го порядка', ...
           'sym_W', K_sym * (T_sym * s + 1), ...
           'tf_val', K_val * (T_val * s_tf + 1) / (1e-4 * s_tf + 1), ...
           'has_ideal_impulse', true), ...
           
    struct('id', 11, 'name', 'Изодромное', ...
           'sym_W', (K_sym * (T_sym * s + 1)) / s, ...
           'tf_val', (K_val * (T_val * s_tf + 1)) / s_tf, ...
           'has_ideal_impulse', false)
};

for i = 1:length(links)
    fprintf('[%d] %s\n', links{i}.id, links{i}.name);
    W_s = links{i}.sym_W;
    
    Y_step = W_s * (1 / s);
    h_analytic = ilaplace(Y_step, s, t);
    
    Y_impulse = W_s * 1;
    w_analytic = ilaplace(Y_impulse, s, t);
    
    links{i}.h_analytic = h_analytic;
    links{i}.w_analytic = w_analytic;
    
    fprintf('W(s) = %s\n', char(W_s));
    fprintf('h(t) = %s\n', char(h_analytic));
    fprintf('w(t) = %s\n\n', char(w_analytic));
end

figure('Name', 'Переходные функции h(t) всех звеньев', 'Color', 'w', 'Position', [50, 50, 1200, 800]);
for i = 1:length(links)
    subplot(4, 3, i);
    hold on; grid on;
    
    if ismember(links{i}.id, [2, 3, 10])
        t_sub = time;
        if links{i}.id == 3
            y_h = K_val * ones(size(t_sub));
        elseif links{i}.id == 2
            y_h = zeros(size(t_sub));
            y_h(1) = K_val / 0.005;
        elseif links{i}.id == 10
            y_h = K_val * ones(size(t_sub));
            y_h(1) = K_val + (K_val * T_val) / 0.005;
        end
        plot(t_sub, y_h, 'b', 'LineWidth', 1.5);
    else
        [y_sim, t_sim] = step(links{i}.tf_val, time);
        plot(t_sim, y_sim, 'b', 'LineWidth', 1.5);
    end
    
    title(sprintf('%d. %s', links{i}.id, links{i}.name), 'FontSize', 9, 'Color', 'k');
    xlabel('t, c', 'FontSize', 8, 'Color', 'k');
    ylabel('h(t)', 'FontSize', 8, 'Color', 'k');
    xlim([0 5]);
    ax = gca;
    ax.Color = [0.93 0.93 0.93];
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.GridColor = [0.6 0.6 0.6];
    ax.GridAlpha = 0.8;
end

figure('Name', 'Импульсные переходные функции w(t) всех звеньев', 'Color', 'w', 'Position', [100, 50, 1200, 800]);
for i = 1:length(links)
    subplot(4, 3, i);
    hold on; grid on;
    
    if links{i}.has_ideal_impulse
        t_sub = time;
        if links{i}.id == 3
            y_w = zeros(size(t_sub));
            y_w(1) = K_val / 0.005;
            plot(t_sub, y_w, 'r', 'LineWidth', 1.5);
        elseif links{i}.id == 2
            y_w = zeros(size(t_sub));
            y_w(1) = K_val / (0.005^2);
            y_w(2) = -K_val / (0.005^2);
            plot(t_sub, y_w, 'r', 'LineWidth', 1.5);
        elseif links{i}.id == 10
            y_w = zeros(size(t_sub));
            y_w(1) = (K_val * T_val) / (0.005^2);
            y_w(2) = -(K_val * T_val) / (0.005^2);
            plot(t_sub, y_w, 'r', 'LineWidth', 1.5);
        end
    else
        [y_sim, t_sim] = impulse(links{i}.tf_val, time);
        plot(t_sim, y_sim, 'r', 'LineWidth', 1.5);
    end
    
    title(sprintf('%d. %s', links{i}.id, links{i}.name), 'FontSize', 9, 'Color', 'k');
    xlabel('t, c', 'FontSize', 8, 'Color', 'k');
    ylabel('w(t)', 'FontSize', 8, 'Color', 'k');
    xlim([0 5]);
    ax = gca;
    ax.Color = [0.93 0.93 0.93];
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.GridColor = [0.6 0.6 0.6];
    ax.GridAlpha = 0.8;
end

figure('Name', 'Влияние коэффициента K', 'Color', 'w');
hold on; grid on;
K_test_vals = [1.0, 2.0, K_val, 4.0, 5.0];
for kv = K_test_vals
    W_k = kv / (T2_val^2 * s_tf^2 + T1_val * s_tf + 1);
    [y_k, t_k] = step(W_k, time);
    plot(t_k, y_k, 'LineWidth', 1.5);
end
title('Колебательное звено: варьирование коэффициента передачи K', 'Color', 'k');
xlabel('Время t, с', 'Color', 'k'); ylabel('h(t)', 'Color', 'k');
legend(arrayfun(@(k) sprintf('K = %.2f', k), K_test_vals, 'UniformOutput', false), 'Location', 'northeast');
ax = gca;
ax.Color = [0.93 0.93 0.93]; ax.XColor = 'k'; ax.YColor = 'k';
ax.GridColor = [0.6 0.6 0.6]; ax.GridAlpha = 0.8;
hold off;

T_factors = [0.2, 0.5, 1.0, 2.0, 5.0];
figure('Name', 'Исследование апериодического звена 1-го порядка при изменении T', 'Color', 'w');
hold on; grid on;

for f = T_factors
    cur_T = T_val * f;
    W_aper = K_val / (cur_T * s_tf + 1);
    [y_ap, t_ap] = step(W_aper, time);
    plot(t_ap, y_ap, 'LineWidth', 1.5);
end
title(sprintf('Апериодическое звено №4: изменение T в 2 и 5 раз (номинал T = %.2f с)', T_val), 'Color', 'k');
xlabel('Время t, с', 'Color', 'k'); ylabel('h(t)', 'Color', 'k');
legend({'T / 5 (0.06 с)', 'T / 2 (0.15 с)', 'Номинал T (0.30 с)', 'T * 2 (0.60 с)', 'T * 5 (1.50 с)'}, 'Location', 'southeast');
ax = gca;
ax.Color = [0.93 0.93 0.93]; ax.XColor = 'k'; ax.YColor = 'k';
ax.GridColor = [0.6 0.6 0.6]; ax.GridAlpha = 0.8;
hold off;

T1_series = [0.8, 0.6, 0.4, 0.2, 0.05];
figure('Name', 'Колебательное звено при уменьшении T1', 'Color', 'w');
hold on; grid on;

for t1 = T1_series
    W_koleb = K_val / (T2_val^2 * s_tf^2 + t1 * s_tf + 1);
    [y_kol, t_kol] = step(W_koleb, 0:0.005:15);
    plot(t_kol, y_kol, 'LineWidth', 1.5);
end
title(sprintf('Колебательное звено №6: уменьшение T_1 при фиксированном T_2 = %.2f с', T2_val), 'Color', 'k');
xlabel('Время t, с', 'Color', 'k'); ylabel('h(t)', 'Color', 'k');
legend(arrayfun(@(t1) sprintf('T_1 = %.2f с (\\xi = %.2f)', t1, t1 / (2 * T2_val)), T1_series, 'UniformOutput', false), 'Location', 'northeast');
ax = gca;
ax.Color = [0.93 0.93 0.93]; ax.XColor = 'k'; ax.YColor = 'k';
ax.GridColor = [0.6 0.6 0.6]; ax.GridAlpha = 0.8;
hold off;

fprintf('Расчеты завершены. Графики построены.\n');
