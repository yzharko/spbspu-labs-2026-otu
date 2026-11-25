clear; clc; close all;

syms s t K;
K0 = 1.6;
T = 0.8;
T1 = 1.0;
T2 = 1.1;
T1_aper = 2.5;
s_tf = tf('s');

links = {
    struct('name', 'Интегрирующее', 'W', K/s, 'W_check', @(k) k/s_tf)
    struct('name', 'Идеальное дифференцирующее', 'W', K*s, 'W_check', @(k) k*s_tf)
    struct('name', 'Усилительное', 'W', K, 'W_check', @(k) tf(k))
    struct('name', 'Апериодическое 1-го порядка', 'W', K/(T*s+1), 'W_check', @(k) k/(T*s_tf+1))
    struct('name', 'Апериодическое 2-го порядка', 'W', K/(T2^2*s^2+T1_aper*s+1), 'W_check', @(k) k/(T2^2*s_tf^2+T1_aper*s_tf+1))
    struct('name', 'Колебательное', 'W', K/(T2^2*s^2+T1*s+1), 'W_check', @(k) k/(T2^2*s_tf^2+T1*s_tf+1))
    struct('name', 'Консервативное', 'W', K/(T*s^2+1), 'W_check', @(k) k/(T*s_tf^2+1))
    struct('name', 'Инерционное интегрирующее', 'W', K/(s*(T*s+1)), 'W_check', @(k) k/(s_tf*(T*s_tf+1)))
    struct('name', 'Реальное дифференцирующее', 'W', K*s/(T*s+1), 'W_check', @(k) k*s_tf/(T*s_tf+1))
    struct('name', 'Форсирующее 1-го порядка', 'W', K*(T*s+1), 'W_check', @(k) k*(T*s_tf+1))
    struct('name', 'Изодромное', 'W', K*(T*s+1)/s, 'W_check', @(k) k*(T*s_tf+1)/s_tf)
};
N = numel(links);
improper = [2 10];

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

for i = 1:N
    Y_step = links{i}.W * U_step;
    links{i}.h = ilaplace(Y_step, s, t);
    Y_impulse = links{i}.W * U_impulse;
    links{i}.w = ilaplace(Y_impulse, s, t);

    fprintf('%d. %s\n', i, links{i}.name);
    disp(['W(s) = ' char(links{i}.W)]);
    disp(['h(t) = ' char(links{i}.h)]);
    disp(['w(t) = ' char(links{i}.w)]);
end

time = 0:0.01:10;
dt = time(2) - time(1);

kinds = {'h', 'Переходные характеристики динамических звеньев'
         'w', 'Импульсные характеристики динамических звеньев'};
for c = 1:size(kinds, 1)
    fld = kinds{c, 1};
    figure;
    tiledlayout(4, 3);
    for i = 1:N
        nexttile;
        grid on;
        f_k = subs(links{i}.(fld), K, K0);
        f_k = subs(f_k, {heaviside(t), dirac(t), dirac(1, t)}, {1, 0, 0});
        y = double(subs(f_k, t, time)) .* ones(size(time));
        plot(time, y);
        title(sprintf('%d. %s', i, links{i}.name));
        xlabel('Время t, с');
        ylabel([fld '(t)']);
    end
    sgtitle(kinds{c, 2});
end

figure;
tiledlayout(4, 3);
for i = 1:N
    nexttile;
    grid on;
    switch i
        case 2
            y = zeros(size(time));
            y(1) = K0/dt;
            plot(time, y);
        case 10
            y = K0*ones(size(time));
            y(1) = K0 + K0*T/dt;
            plot(time, y);
        otherwise
            W_check = links{i}.W_check(K0);
            [y, t_out] = step(W_check, time);
            plot(t_out, y);
    end
    title(sprintf('%d. %s', i, links{i}.name));
    xlabel('Время t, с');
    ylabel('Выходной сигнал');
end
sgtitle('Проверка переходной характеристики');

figure;
tiledlayout('flow');
for i = setdiff(1:N, improper)
    nexttile;
    grid on;
    W_check = links{i}.W_check(K0);
    [y, t_out] = impulse(W_check, time);
    plot(t_out, y);
    title(sprintf('%d. %s', i, links{i}.name));
    xlabel('Время t, с');
    ylabel('Выходной сигнал');
end
sgtitle('Проверка импульсной характеристики');

K_values = [0.5 1 1.6 3];
figure;
tiledlayout(1, 2);
for i = [4 6]
    nexttile;
    hold on;
    grid on;
    for k = K_values
        W_check = links{i}.W_check(k);
        plot(time, step(W_check, time), 'DisplayName', ['K = ', num2str(k)]);
    end
    title(sprintf('Звено %d: h(t) при изменении K', i));
    xlabel('Время t, с');
    ylabel('h(t)');
    legend('show', 'Location', 'southeast');
end

time_T = 0:0.01:25;
T_values = T*[1/5 1/2 1 2 5];
figure;
hold on;
grid on;
fprintf('Звено 4: влияние T при K = %g\n', K0);
for Tk = T_values
    W_check = tf(K0, [Tk 1]);
    plot(time_T, step(W_check, time_T), 'DisplayName', ['T = ', num2str(Tk)]);
    info = stepinfo(W_check, 'SettlingTimeThreshold', 0.05);
    fprintf('T = %-5g tп(5%%) = %5.2f c   tп/T = %.2f\n', Tk, info.SettlingTime, info.SettlingTime/Tk);
end
title('Звено 4: h(t) при изменении T');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'southeast');

time_T1 = 0:0.05:100;
T1_values = [1 0.5 0.25 0.1];
figure;
hold on;
grid on;
fprintf('Звено 6: влияние T1 при T2 = %g, K = %g\n', T2, K0);
for T1k = T1_values
    W_check = tf(K0, [T2^2 T1k 1]);
    plot(time_T1, step(W_check, time_T1), 'DisplayName', ['T1 = ', num2str(T1k)]);
    info = stepinfo(W_check, 'SettlingTimeThreshold', 0.05);
    fprintf('T1 = %-5g zeta = %.3f   перерегулирование = %5.2f %%   tп(5%%) = %6.2f c\n', T1k, T1k/(2*T2), info.Overshoot, info.SettlingTime);
end
title('Звено 6: h(t) при уменьшении T1');
xlabel('Время t, с');
ylabel('h(t)');
legend('show', 'Location', 'southeast');

