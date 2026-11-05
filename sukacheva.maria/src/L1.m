clear;
clc;
close all;

syms s t K;

K_val = 1.6;
T = 0.8;
T1 = 1.0;
T1_2_aperiodic = 3.0;
T2 = 1.1;

W_integrial_ideal = K / s;
W_diff_ideal = K * s;
W_noninertial_amplifying = K;
W_1_aperiodic = K / (T * s + 1);
W_2_aperiodic = K / (T2^2 * s^2 + T1_2_aperiodic * s + 1);
W_oscillatory = K / (T2^2 * s^2 + T1 * s + 1);
W_conservative = K / (T * s^2 + 1);
W_integrial_inertial = K / (s * (T * s + 1));
W_diff_real = (K * s) / (T * s + 1);
W_1_forcing = K * (T * s + 1);
W_isodromic = K * (T * s + 1) / s;

U_step = laplace(heaviside(t), t, s);

Y_step_integrial_ideal = W_integrial_ideal * U_step;
Y_step_diff_ideal = W_diff_ideal * U_step;
Y_step_noninertial_amplifying = W_noninertial_amplifying * U_step;
Y_step_1_aperiodic = W_1_aperiodic * U_step;
Y_step_2_aperiodic = W_2_aperiodic * U_step;
Y_step_oscillatory = W_oscillatory * U_step;
Y_step_conservative = W_conservative * U_step;
Y_step_integrial_inertial = W_integrial_inertial * U_step;
Y_step_diff_real = W_diff_real * U_step;
Y_step_1_forcing = W_1_forcing * U_step;
Y_step_isodromic = W_isodromic * U_step;

h_integrial_ideal = ilaplace(Y_step_integrial_ideal, s, t);
h_diff_ideal = ilaplace(Y_step_diff_ideal, s, t);
h_noninertial_amplifying = ilaplace(Y_step_noninertial_amplifying, s, t);
h_1_aperiodic = ilaplace(Y_step_1_aperiodic, s, t);
h_2_aperiodic = ilaplace(Y_step_2_aperiodic, s, t);
h_oscillatory = ilaplace(Y_step_oscillatory, s, t);
h_conservative = ilaplace(Y_step_conservative, s, t);
h_integrial_inertial = ilaplace(Y_step_integrial_inertial, s, t);
h_diff_real = ilaplace(Y_step_diff_real, s, t);
h_1_forcing = ilaplace(Y_step_1_forcing, s, t);
h_isodromic = ilaplace(Y_step_isodromic, s, t);

U_impulse = laplace(dirac(t), t, s);

Y_impulse_integrial_ideal = W_integrial_ideal * U_impulse;
Y_impulse_diff_ideal = W_diff_ideal * U_impulse;
Y_impulse_noninertial_amplifying = W_noninertial_amplifying * U_impulse;
Y_impulse_1_aperiodic = W_1_aperiodic * U_impulse;
Y_impulse_2_aperiodic = W_2_aperiodic * U_impulse;
Y_impulse_oscillatory = W_oscillatory * U_impulse;
Y_impulse_conservative = W_conservative * U_impulse;
Y_impulse_integrial_inertial = W_integrial_inertial * U_impulse;
Y_impulse_diff_real = W_diff_real * U_impulse;
Y_impulse_1_forcing = W_1_forcing * U_impulse;
Y_impulse_isodromic = W_isodromic * U_impulse;

w_integrial_ideal = ilaplace(Y_impulse_integrial_ideal, s, t);
w_diff_ideal = ilaplace(Y_impulse_diff_ideal, s, t);
w_noninertial_amplifying = ilaplace(Y_impulse_noninertial_amplifying, s, t);
w_1_aperiodic = ilaplace(Y_impulse_1_aperiodic, s, t);
w_2_aperiodic = ilaplace(Y_impulse_2_aperiodic, s, t);
w_oscillatory = ilaplace(Y_impulse_oscillatory, s, t);
w_conservative = ilaplace(Y_impulse_conservative, s, t);
w_integrial_inertial = ilaplace(Y_impulse_integrial_inertial, s, t);
w_diff_real = ilaplace(Y_impulse_diff_real, s, t);
w_1_forcing = ilaplace(Y_impulse_1_forcing, s, t);
w_isodromic = ilaplace(Y_impulse_isodromic, s, t);

names = {
    '1. Интегрирующее'
    '2. Дифференцирующее (идеальное)'
    '3. Усилительное (безынерционное)'
    '4. Апериодическое 1-го порядка'
    '5. Апериодическое 2-го порядка'
    '6. Колебательное'
    '7. Консервативное'
    '8. Инерционное интегрирующее'
    '9. Реальное дифференцирующее'
    '10. Форсирующее 1-го порядка'
    '11. Изодромное'
};

h_all = {
    h_integrial_ideal
    h_diff_ideal
    h_noninertial_amplifying
    h_1_aperiodic
    h_2_aperiodic
    h_oscillatory
    h_conservative
    h_integrial_inertial
    h_diff_real
    h_1_forcing
    h_isodromic
};

w_all = {
    w_integrial_ideal
    w_diff_ideal
    w_noninertial_amplifying
    w_1_aperiodic
    w_2_aperiodic
    w_oscillatory
    w_conservative
    w_integrial_inertial
    w_diff_real
    w_1_forcing
    w_isodromic
};

time = 0:0.01:10;
time_osc = 0:0.01:40;

plot_all(h_all, names, K, K_val, t, time, 'Переходные характеристики h(t)', 'h(t)');
plot_all(w_all, names, K, K_val, t, time, 'Импульсные характеристики w(t)', 'w(t)');

s_tf = tf('s');

Wc_all = {
    K_val / s_tf
    K_val * s_tf
    tf(K_val)
    K_val / (T * s_tf + 1)
    K_val / (T2^2 * s_tf^2 + T1_2_aperiodic * s_tf + 1)
    K_val / (T2^2 * s_tf^2 + T1 * s_tf + 1)
    K_val / (T * s_tf^2 + 1)
    K_val / (s_tf * (T * s_tf + 1))
    (K_val * s_tf) / (T * s_tf + 1)
    K_val * (T * s_tf + 1)
    (K_val * (T * s_tf + 1)) / s_tf
};

check_all(Wc_all, names, time, 'step', 'Проверка переходных характеристик (Toolbox)');
check_all(Wc_all, names, time, 'impulse', 'Проверка импульсных характеристик (Toolbox)');

figure;
hold on;
grid on;

for m = [1/5 1/2 1 2 5]
    Wt = K_val / (T * m * s_tf + 1);
    plot(time, step(Wt, time), 'DisplayName', ['T = ', num2str(T * m)]);
end

title('Апериодическое звено: влияние T');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'southeast');
hold off;

figure;
hold on;
grid on;

for T1k = [1.9 1.5 1.0 0.7 0.4 0.2 0.05]
    Wt = K_val / (T2^2 * s_tf^2 + T1k * s_tf + 1);
    plot(time_osc, step(Wt, time_osc), 'DisplayName', ['T1 = ', num2str(T1k)]);
end

title('Колебательное звено: влияние T1 при T2 = 1.1');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'northeast');
hold off;

function plot_all(expr_list, names, K, K_val, t, time, ttl, yl)

figure('Position', [100 100 1300 900]);

for i = 1:numel(expr_list)

    subplot(4, 3, i);
    hold on;
    grid on;

    expr = subs(expr_list{i}, K, K_val);
    expr = subs(expr, dirac(t, 1), 0);
    expr = subs(expr, dirac(t), 0);

    y = double(subs(expr, t, time));
    y = real(y(:));

    plot(time, y, 'LineWidth', 1.5);

    title(names{i}, 'FontSize', 9);
    xlabel('Время t, с');
    ylabel(yl);

    hold off;
end

subplot(4, 3, 12);
axis off;
text(0.1, 0.6, ['K = ', num2str(K_val)], 'FontSize', 13);

sgtitle(ttl);

end

function check_all(Wc_all, names, time, kind, ttl)

figure('Position', [100 100 1300 900]);

for i = 1:numel(Wc_all)

    subplot(4, 3, i);
    hold on;
    grid on;

    try
        if strcmp(kind, 'step')
            y = step(Wc_all{i}, time);
        else
            y = impulse(Wc_all{i}, time);
        end

        plot(time, y(:), 'LineWidth', 1.5);
        title(names{i}, 'FontSize', 9);

    catch
        title([names{i}, ' (не вычисляется)'], 'FontSize', 9);
    end

    xlabel('Время t, с');
    ylabel('Выходной сигнал');

    hold off;
end

sgtitle(ttl);

end
