clear;
clc;
close all;

K0  = 2.0;
T0  = 0.5;
T10 = 0.5;
T20 = 0.4;
T1a = 2*T20;
time = 0:0.01:10;
syms s t K
s_tf = tf('s');

names = {'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
    'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', ...
    'Колебательное', 'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее', 'Изодромное'};

W = {K/s, K*s, K, K/(T0*s+1), K/(T20^2*s^2+T1a*s+1), ...
    K/(T20^2*s^2+T10*s+1), K/(T0*s^2+1), K/(s*(T0*s+1)), ...
    K*s/(T0*s+1), K*(T0*s+1), K*(T0*s+1)/s};

U_step = laplace(heaviside(t), t, s);
U_imp  = laplace(dirac(t), t, s);
h = cell(1, 11); w = cell(1, 11);
for i = 1:11
    h{i} = simplify(ilaplace(W{i}*U_step, s, t));
    w{i} = simplify(ilaplace(W{i}*U_imp, s, t));
    fprintf('%d. %s\n  W(s) = %s\n  h(t) = %s\n  w(t) = %s\n\n', ...
        i, names{i}, char(W{i}), char(h{i}), char(w{i}));
end

for k = 1:2
    figure;
    for i = 1:11
        subplot(3, 4, i); hold on; grid on;
        if k == 1
            expr = subs(h{i}, K, K0);
        else
            expr = subs(w{i}, K, K0);
        end
        expr = subs(expr, {heaviside(t), dirac(t), dirac(1, t)}, {1, 0, 0});
        y = double(subs(expr, t, time));
        plot(time, y, 'LineWidth', 1.4);
        if k == 1
            ylabel('h(t)');
        else
            ylabel('w(t)');
        end
        title(sprintf('%d. %s', i, names{i}));
        xlabel('t, с');
    end
end

W_tf = {K0/s_tf, K0*s_tf, tf(K0, 1), K0/(T0*s_tf+1), ...
    K0/(T20^2*s_tf^2+T1a*s_tf+1), K0/(T20^2*s_tf^2+T10*s_tf+1), ...
    K0/(T0*s_tf^2+1), K0/(s_tf*(T0*s_tf+1)), K0*s_tf/(T0*s_tf+1), ...
    K0*(T0*s_tf+1), K0*(T0*s_tf+1)/s_tf};
for i = 1:11
    h_num = double(subs(subs(h{i}, K, K0), t, time));
    try
        y = step(W_tf{i}, time);
        fprintf('Звено %2d: max|h - h_ЧС| = %.2e\n', i, max(abs(h_num(:)-y(:))));
    catch
        fprintf('Звено %2d: проверка step неприменима (несобственная ПФ)\n', i);
    end
end

K_models = {@(Kv) tf(Kv, [1 0]), @(Kv) tf(Kv, [T0 1]), @(Kv) tf(Kv, [T20^2 T10 1])};
K_titles = {'Интегрирующее звено', 'Апериодическое звено 1-го порядка', 'Колебательное звено'};
for j = 1:3
    figure; hold on; grid on;
    for Kv = [0.5 1 2 4]
        plot(time, step(K_models{j}(Kv), time), 'DisplayName', sprintf('K = %g', Kv));
    end
    title(K_titles{j});
    legend('show'); xlabel('t, с'); ylabel('h(t)');
end

fprintf('\nВлияние T:\n');
T_values = T0*[1/5 1/2 1 2 5];
figure; hold on; grid on;
for Tv = T_values
    S = stepinfo(tf(K0, [Tv 1]));
    fprintf('T = %.2f с: t_рег = %.3f с\n', Tv, S.SettlingTime);
    plot(time, step(tf(K0, [Tv 1]), time), 'DisplayName', sprintf('T = %.2f', Tv));
end
legend('show'); xlabel('t, с'); ylabel('h(t)');

fprintf('\nВлияние T1 (T2 = %.2f):\n', T20);
T1_values = [0.5 0.4 0.3 0.2 0.1 0.05];
figure; hold on; grid on;
for T1v = T1_values
    S = stepinfo(tf(K0, [T20^2 T1v 1]));
    fprintf('T1 = %.2f с: перерегулирование = %.2f %%, t_зат = %.3f с\n', ...
        T1v, S.Overshoot, S.SettlingTime);
    plot(time, step(tf(K0, [T20^2 T1v 1]), time), 'DisplayName', sprintf('T1 = %.2f', T1v));
end
legend('show'); xlabel('t, с'); ylabel('h(t)');
