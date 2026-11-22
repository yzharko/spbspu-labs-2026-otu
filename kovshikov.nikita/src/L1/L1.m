syms s t K;
K_val = 3.0;
T = 1.0;
T1 = 1.2;
T2 = 1.5;
T1_aperiodic = 3.6;

K_values = 1:10;
time = 0:0.01:10;

W_names = {'Интегрирующее', 'Идеальное дифференцирующее', ...
    'Усилительное', 'Апериодическое 1-го порядка', ...
    'Апериодическое 2-го порядка', 'Колебательное', ...
    'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее 1-го порядка', 'Изодромное'};

W_list = {K/s, K*s, K, K/(T*s+1), ...
    K/(T2^2*s^2+T1_aperiodic*s+1), K/(T2^2*s^2+T1*s+1), ...
    K/(T*s^2+1), K/(s*(T*s+1)), K*s/(T*s+1), ...
    K*(T*s+1), K*(T*s+1)/s};

s_tf = tf('s');
W_check_list = {1/s_tf, s_tf, tf(1), 1/(T*s_tf+1), ...
    1/(T2^2*s_tf^2+T1_aperiodic*s_tf+1), ...
    1/(T2^2*s_tf^2+T1*s_tf+1), 1/(T*s_tf^2+1), ...
    1/(s_tf*(T*s_tf+1)), s_tf/(T*s_tf+1), ...
    T*s_tf+1, (T*s_tf+1)/s_tf};

for i = 1:numel(W_list)
    W = W_list{i};
    U_step = laplace(heaviside(t), t, s);
    Y_step = W * U_step;
    h = ilaplace(Y_step, s, t);
    U_impulse = laplace(dirac(t), t, s);
    Y_impulse = W * U_impulse;
    w = ilaplace(Y_impulse, s, t);

    %% 1. Построение переходных характеристик
    figure;
    hold on;
    grid on;
    for k = K_values
        h_k = subs(h, K, k);
        y = double(subs(h_k, t, time));
        y = y + zeros(size(time));
        plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
    end
    title([W_names{i}, ': переходные хар-ки']);
    xlabel('Время t, с');
    ylabel('h(t)');
    legend('show', 'Location', 'best');
    hold off;

    %% 2. Построение импульсных характеристик
    figure;
    hold on;
    grid on;
    for k = K_values
        w_k = subs(w, K, k);
        y = double(subs(w_k, t, time));
        y = y + zeros(size(time));
        plot(time, y, 'DisplayName', ['K = ', num2str(k)]);
    end
    title([W_names{i}, ': импульсные хар-ки']);
    xlabel('Время t, с');
    ylabel('w(t)');
    legend('show', 'Location', 'northwest');
    hold off;

    %% 3. Проверка с помощью Control System Toolbox
    if isproper(W_check_list{i})
        figure;
        hold on;
        grid on;
        for k = K_values
            W_check = k * W_check_list{i};
            [y, tt] = step(W_check, time);
            plot(tt, y, 'DisplayName', ['K = ', num2str(k)]);
        end
        title([W_names{i}, ': проверка step']);
        xlabel('Время t, с');
        ylabel('h(t)');
        legend('show', 'Location', 'northwest');
        hold off;

        figure;
        hold on;
        grid on;
        for k = K_values
            W_check = k * W_check_list{i};
            [y, tt] = impulse(W_check, time);
            plot(tt, y, 'DisplayName', ['K = ', num2str(k)]);
        end
        title([W_names{i}, ': проверка impulse']);
        xlabel('Время t, с');
        ylabel('w(t)');
        legend('show', 'Location', 'northwest');
        hold off;
    end
end

%% 4. изменение T в 2 и 5 раз
T_values = T * [1/5, 1/2, 1, 2, 5];
time_T = 0:0.01:30;
figure;
hold on;
grid on;
for T_current = T_values
    W_current = K_val / (T_current*s_tf+1);
    [y, tt] = step(W_current, time_T);
    plot(tt, y, 'DisplayName', ['T = ', num2str(T_current)]);
end
title('Влияние T на апериодическое звено при K = 3');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'northwest');
hold off;

%% 5. уменьшение T1 при фиксированном T2
T1_values = T1 * [1, 1/2, 1/5];
time_long = 0:0.01:100;
figure;
hold on;
grid on;
for T1_current = T1_values
    W_current = K_val / (T2^2*s_tf^2+T1_current*s_tf+1);
    [y, tt] = step(W_current, time_long);
    plot(tt, y, 'DisplayName', ['T1 = ', num2str(T1_current)]);
end
title('Влияние T1 на колебательное звено: K = 3, T2 = 1.5');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'northwest');
hold off;

