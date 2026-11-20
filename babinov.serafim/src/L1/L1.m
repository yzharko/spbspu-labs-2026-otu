clear; clc; close all;

K = 2.1;
T = 0.9;
T1 = 1.0;
T2 = 1.1;
T1a = 2*T2;

time = 0:0.01:20;
names = { ...
    'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
    'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', ...
    'Колебательное', 'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее', 'Изодромное'};

syms s t k real

Wsym = { ...
    k/s, ...
    k*s, ...
    k, ...
    k/(T*s + 1), ...
    k/(T2^2*s^2 + T1a*s + 1), ...
    k/(T2^2*s^2 + T1*s + 1), ...
    k/(T*s^2 + 1), ...
    k/(s*(T*s + 1)), ...
    k*s/(T*s + 1), ...
    k*(T*s + 1), ...
    k*(T*s + 1)/s};

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

for i = 1:numel(Wsym)
    W = subs(Wsym{i}, k, K);
    h = simplify(ilaplace(W*U_step, s, t));
    w = simplify(ilaplace(W*U_impulse, s, t));

    fprintf('\n%d. %s\n', i, names{i});
    fprintf('W(s) = %s\n', char(W));
    fprintf('h(t) = %s\n', char(h));
    fprintf('w(t) = %s\n', char(w));
end

s_tf = tf('s');
W = { ...
    K/s_tf, ...
    K*s_tf, ...
    tf(K), ...
    K/(T*s_tf + 1), ...
    K/(T2^2*s_tf^2 + T1a*s_tf + 1), ...
    K/(T2^2*s_tf^2 + T1*s_tf + 1), ...
    K/(T*s_tf^2 + 1), ...
    K/(s_tf*(T*s_tf + 1)), ...
    K*s_tf/(T*s_tf + 1), ...
    K*(T*s_tf + 1), ...
    K*(T*s_tf + 1)/s_tf};

figure('Name', 'Переходные характеристики');
for i = 1:numel(W)
    subplot(3, 4, i); hold on; grid on;
    if i == 2
        y = zeros(size(time));
        plot(time, y);
        stem(0, K, 'b', 'filled');
    elseif i == 10
        y = K*ones(size(time));
        plot(time, y);
        stem(0, K*T, 'b', 'filled');
    else
        [y, t_out] = step(W{i}, time);
        plot(t_out, y);
    end
    title(sprintf('%d. %s', i, names{i}), 'FontSize', 9);
    xlabel('t, с'); ylabel('h(t)');
    hold off;
end
sgtitle('Переходные характеристики типовых звеньев');

figure('Name', 'Импульсные характеристики');
for i = 1:numel(W)
    subplot(3, 4, i); hold on; grid on;
    if i == 2
        y = zeros(size(time));
        plot(time, y);
        stem([-time(2) time(2)], [K -K], 'b', 'filled');
    elseif i == 3
        y = zeros(size(time));
        plot(time, y);
        stem(0, K, 'b', 'filled');
    elseif i == 10
        y = zeros(size(time));
        plot(time, y);
        stem([-time(2) 0 time(2)], [K*T K -K*T], 'b', 'filled');
    elseif i == 9
        y = -(K/T^2)*exp(-time/T);
        plot(time, y);
    elseif i == 11
        y = K*ones(size(time));
        plot(time, y);
    else
        [y, t_out] = impulse(W{i}, time);
        plot(t_out, y);
    end
    title(sprintf('%d. %s', i, names{i}), 'FontSize', 9);
    xlabel('t, с'); ylabel('w(t)');
    hold off;
end
sgtitle('Импульсные характеристики типовых звеньев');

figure('Name', 'Влияние K'); hold on; grid on;
for Ki = 1:10
    [y, t_out] = step(Ki/(T2^2*s_tf^2 + T1*s_tf + 1), time);
    plot(t_out, y);
end
title('Колебательное звено при K = 1...10');
xlabel('t, с'); ylabel('h(t)');
legend(compose('K = %d', 1:10), 'Location', 'best');
hold off;

T_values = T*[0.2 0.5 1 2 5];
figure('Name', 'Влияние T'); hold on; grid on;
for Ti = T_values
    [y, t_out] = step(K/(Ti*s_tf + 1), time);
    plot(t_out, y);
end
title('Апериодическое звено 1-го порядка при разных T');
xlabel('t, с'); ylabel('h(t)');
legend(compose('T = %.2g', T_values), 'Location', 'best');
hold off;

T1_values = [T1 0.75 0.5 0.25 0.1];
figure('Name', 'Влияние T1'); hold on; grid on;
for T1i = T1_values
    [y, t_out] = step(K/(T2^2*s_tf^2 + T1i*s_tf + 1), time);
    plot(t_out, y);
end
title('Колебательное звено при уменьшении T1');
xlabel('t, с'); ylabel('h(t)');
legend(compose('T1 = %.2g', T1_values), 'Location', 'best');
hold off;
