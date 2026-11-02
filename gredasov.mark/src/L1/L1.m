clear;
clc;
close all;

params.K = 3;
params.T = 1;
params.T1 = 1.2;
params.T2 = 1.5;
params.T1a = 3.5;

eps = 0.001;
dt = 0.05;
time = 0:dt:10;

syms s t K T T1 T2 T1a real

s_tf = tf('s');

links = [
    struct('n', 1, ...
           'name', 'Интегрирующее (идеальное)', ...
           'W', K/s, ...
           'W_check', params.K/s_tf)

    struct('n', 2, ...
           'name', 'Дифференцирующее (идеальное)', ...
           'W', K*s, ...
           'W_check', params.K*s_tf/(eps*s_tf + 1))

    struct('n', 3, ...
           'name', 'Усилительное (безынерционное)', ...
           'W', K, ...
           'W_check', params.K)

    struct('n', 4, ...
           'name', 'Апериодическое 1-го порядка', ...
           'W', K/(T*s + 1), ...
           'W_check', params.K/(params.T*s_tf + 1))

    struct('n', 5, ...
           'name', 'Апериодическое 2-го порядка', ...
           'W', K/(T2^2*s^2 + T1a*s + 1), ...
           'W_check', params.K/(params.T2^2*s_tf^2 + params.T1a*s_tf + 1))

    struct('n', 6, ...
           'name', 'Колебательное', ...
           'W', K/(T2^2*s^2 + T1*s + 1), ...
           'W_check', params.K/(params.T2^2*s_tf^2 + params.T1*s_tf + 1))

    struct('n', 7, ...
           'name', 'Консервативное', ...
           'W', K/(T*s^2 + 1), ...
           'W_check', params.K/(params.T*s_tf^2 + 1))

    struct('n', 8, ...
           'name', 'Инерционное интегрирующее', ...
           'W', K/(s*(T*s + 1)), ...
           'W_check', params.K/(s_tf*(params.T*s_tf + 1)))

    struct('n', 9, ...
           'name', 'Реальное дифференцирующее', ...
           'W', K*s/(T*s + 1), ...
           'W_check', params.K*s_tf/(params.T*s_tf + 1))

    struct('n', 10, ...
           'name', 'Форсирующее 1-го порядка', ...
           'W', K*(T*s + 1), ...
           'W_check', params.K*(params.T*s_tf + 1)/(eps*s_tf + 1))

    struct('n', 11, ...
           'name', 'Изодромное', ...
           'W', K*(T*s + 1)/s, ...
           'W_check', params.K*(params.T*s_tf + 1)/s_tf)
];


U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

fprintf('Единичное ступенчатое воздействие: U_step = %s\n', char(U_step))
fprintf('Импульсное воздействие: U_impulse = %s\n\n', char(U_impulse))

for i = 1:length(links)
    links(i).W_num = subs(links(i).W, ...
                          [K, T, T1, T2, T1a], ...
                          [params.K, params.T, params.T1, params.T2, params.T1a]);

    links(i).Y_step = simplify(links(i).W_num * U_step);
    links(i).h = simplify(ilaplace(links(i).Y_step, s, t));

    links(i).Y_impulse = simplify(links(i).W_num * U_impulse);
    links(i).w = simplify(ilaplace(links(i).Y_impulse, s, t));

    fprintf('Звено №%d: %s\n', links(i).n, links(i).name);
    fprintf('\tW(s) = %s\n', char(links(i).W_num));
    fprintf('\tY_step(s) = %s\n', char(links(i).Y_step));
    fprintf('\th(t) = %s\n', char(links(i).h));
    fprintf('\tY_impulse(s) = %s\n', char(links(i).Y_impulse));
    fprintf('\tw(t) = %s\n\n', char(links(i).w));
end

figure('Name', 'Переходные характеристики', 'Color', 'w');
tiledlayout(3, 4);

for i = 1:length(links)
    nexttile;

    if i == 2
        y_analytical = zeros(size(time));
        y_analytical(1) = params.K/dt;
    elseif i == 3
        y_analytical = params.K*ones(size(time));
    elseif i == 10
        y_analytical = params.K*ones(size(time));
        y_analytical(1) = params.K + params.K*params.T/dt;
    else
        h_fun = matlabFunction(links(i).h, 'Vars', t);
        y_analytical = h_fun(time);
    end

    plot(time, y_analytical, 'b', 'LineWidth', 1.5);
    hold on;

    if i == 3
        t_check = time;
        y_check = params.K * ones(size(time));
    else
        [y_check, t_check] = step(links(i).W_check, time);
    end

    plot(t_check, y_check, 'r--', 'LineWidth', 1.5);

    grid on;
    xlabel('t, с');
    ylabel('h(t)');

    title(sprintf('№%d. %s', links(i).n, links(i).name), 'Interpreter', 'none');

    hold off;
end
sgtitle('Сравнение переходных характеристик');

figure('Name', 'Импульсные характеристики', 'Color', 'w');
tiledlayout(3, 4);

for i = 1:length(links)
    nexttile;

    if i == 1
        y_analytical = params.K*ones(size(time));
    elseif i == 2
        y_analytical = zeros(size(time));
        y_analytical(1) = params.K/dt^2;
        y_analytical(2) = -params.K/dt^2;
    elseif i == 3
        y_analytical = zeros(size(time));
        y_analytical(1) = params.K/dt;
    elseif i == 9
        y_analytical = -params.K/params.T^2 * exp(-time/params.T);
        y_analytical(1) = y_analytical(1) + params.K/(params.T*dt);
    elseif i == 10
        y_analytical = zeros(size(time));
        y_analytical(1) = params.K/dt + params.K*params.T/dt^2;
        y_analytical(2) = -params.K*params.T/dt^2;
    elseif i == 11
        y_analytical = params.K*ones(size(time));
        y_analytical(1) = y_analytical(1) + params.K*params.T/dt;
    else
        w_fun = matlabFunction(links(i).w, 'Vars', t);
        y_analytical = w_fun(time);
    end

    plot(time, y_analytical, 'b', 'LineWidth', 1.5);
    hold on;

    if i == 3
        t_check = time;
        y_check = zeros(size(time));
        y_check(1) = params.K/dt;
    else
        [y_check, t_check] = impulse(links(i).W_check, time);
    end
    
    plot(t_check, y_check, 'r--', 'LineWidth', 1.5);

    grid on;
    xlabel('t, с');
    ylabel('w(t)');
    ylim([-3.5 3.5])

    title(sprintf('№%d. %s', links(i).n, links(i).name), 'Interpreter', 'none');

    hold off;
end

sgtitle('Сравнение импульсных характеристик');

K_values = params.K * [0.5 1 1.5 2 3];
time_T = 0:dt:40;

figure('Name', 'Влияние коэффициента K', 'Color', 'w');
hold on;
grid on;

for i = 1:length(K_values)
    K_current = K_values(i);
    W_current = K_current / (params.T*s_tf + 1);
    [y_current, t_current] = step(W_current, time_T);
    plot(t_current, squeeze(y_current), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('K = %.1f', K_current));
end

xlabel('t, с');
ylabel('h(t)');
title('Влияние коэффициента K на переходный процесс');
legend('Location', 'best');
hold off;

T_values = params.T * [0.2 0.5 1 2 5];

figure('Name', 'Влияние постоянной времени T', 'Color', 'w');
hold on;
grid on;

for i = 1:length(T_values)
    T_current = T_values(i);
    W_current = params.K / (T_current*s_tf + 1);
    [y_current, t_current] = step(W_current, time_T);
    plot(t_current, squeeze(y_current), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('T = %.1f', T_current));
end

xlabel('t, с');
ylabel('h(t)');
title('Влияние постоянной времени T на переходный процесс');
legend('Location', 'best');
hold off;

T1_values = params.T1 * [1 0.5 0.25 0.125];

time_T1 = 0:dt:100;

figure('Name', 'Влияние коэффициента T1', 'Color', 'w');
hold on;
grid on;

for i = 1:length(T1_values)
    T1_current = T1_values(i);
    W_current = params.K / (params.T2^2*s_tf^2 + T1_current*s_tf + 1);
    [y_current, t_current] = step(W_current, time_T1);
    plot(t_current, squeeze(y_current), 'LineWidth', 1.5, ...
         'DisplayName', sprintf('T1 = %.2f', T1_current));
end

xlabel('t, с');
ylabel('h(t)');
title('Влияние коэффициента T1 на переходный процесс');
legend('Location', 'best');
hold off;
