clear; clc; close all;

K_val = 4.0;
T = 1.5;
T1 = 1.6;
T2 = 2.0;
T1a = 2*T2;

time = 0:0.01:40;
d = time(2);
names = { ...
    'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
    'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', ...
    'Колебательное', 'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее', 'Изодромное'};

syms s t K

W_sym = { ...
    K/s, K*s, K, K/(T*s + 1), ...
    K/(T2^2*s^2 + T1a*s + 1), K/(T2^2*s^2 + T1*s + 1), ...
    K/(T*s^2 + 1), K/(s*(T*s + 1)), K*s/(T*s + 1), ...
    K*(T*s + 1), K*(T*s + 1)/s};

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

hw = cell(2, numel(W_sym));
for i = 1:numel(W_sym)
    Y_step = W_sym{i}*U_step;
    hw{1, i} = simplify(ilaplace(Y_step, s, t));
    Y_impulse = W_sym{i}*U_impulse;
    hw{2, i} = simplify(ilaplace(Y_impulse, s, t));
    fprintf('\n%d. %s\nW(s) = %s\nh(t) = %s\nw(t) = %s\n', ...
        i, names{i}, char(W_sym{i}), char(hw{1, i}), char(hw{2, i}));
end

s_tf = tf('s');
W_tf = { ...
    K_val/s_tf, K_val*s_tf, tf(K_val), K_val/(T*s_tf + 1), ...
    K_val/(T2^2*s_tf^2 + T1a*s_tf + 1), K_val/(T2^2*s_tf^2 + T1*s_tf + 1), ...
    K_val/(T*s_tf^2 + 1), K_val/(s_tf*(T*s_tf + 1)), K_val*s_tf/(T*s_tf + 1), ...
    K_val*(T*s_tf + 1), K_val*(T*s_tf + 1)/s_tf};

z = zeros(size(time));
c = K_val*ones(size(time));
spec = cell(2, numel(W_tf));
spec{1, 2} = {z, 0, K_val};
spec{1, 10} = {c, 0, K_val*T};
spec{2, 2} = {z, [-d d], [K_val -K_val]};
spec{2, 3} = {z, 0, K_val};
spec{2, 9} = {-(K_val/T^2)*exp(-time/T), 0, K_val/T};
spec{2, 10} = {z, [-d 0 d], [K_val*T K_val -K_val*T]};
spec{2, 11} = {c, 0, K_val*T};

fn = {@step, @impulse};
lab = {'h(t)', 'w(t)'};
kind = {'Переходные', 'Импульсные'};

for p = 1:2
    ttl = [kind{p} ' характеристики типовых звеньев'];
    figure('Name', ttl);
    for i = 1:numel(W_tf)
        subplot(3, 4, i); hold on; grid on;
        if isempty(spec{p, i})
            [y, t_out] = fn{p}(W_tf{i}, time);
            plot(t_out, y);
        else
            plot(time, spec{p, i}{1});
            stem(spec{p, i}{2}, spec{p, i}{3}, 'b', 'filled');
        end
        title(sprintf('%d. %s', i, names{i}), 'FontSize', 9);
        xlabel('t, с'); ylabel(lab{p});
    end
    sgtitle(ttl);
end

K_values = 1:10;
time_k = 0:0.01:10;

ttl = 'Переходные характеристики интегрирующего звена';
figure('Name', ttl); hold on; grid on;
for k = K_values
    y = double(subs(subs(hw{1, 1}, K, k), t, time_k)).*ones(size(time_k));
    plot(time_k, y, 'DisplayName', ['K = ', num2str(k)]);
end
title(ttl);
xlabel('t, с'); ylabel('h(t)');
legend('show', 'Location', 'northwest');

T_values = T*[0.2 0.5 1 2 5];
T1_values = [T1 0.75 0.5 0.25 0.1 0];
sweeps = { ...
    'Переходные характеристики апериодического звена при разных T', ...
    arrayfun(@(x) K_val/(x*s_tf + 1), T_values, 'UniformOutput', false), ...
    compose('T = %.2g', T_values); ...
    'Переходные характеристики колебательного звена при уменьшении T1', ...
    arrayfun(@(x) K_val/(T2^2*s_tf^2 + x*s_tf + 1), T1_values, 'UniformOutput', false), ...
    compose('T1 = %.2g', T1_values)};

for j = 1:2
    figure('Name', sweeps{j, 1}); hold on; grid on;
    for q = 1:numel(sweeps{j, 3})
        [y, t_out] = step(sweeps{j, 2}{q}, time);
        plot(t_out, y);
        info = stepinfo(sweeps{j, 2}{q});
        if j == 1
            fprintf('%s: время нарастания %.3f с, время установления %.3f с\n', ...
                sweeps{j, 3}{q}, info.RiseTime, info.SettlingTime);
        else
            fprintf('%s: перерегулирование %.2f %%, время установления %.3f с\n', ...
                sweeps{j, 3}{q}, info.Overshoot, info.SettlingTime);
        end
    end
    title(sweeps{j, 1});
    xlabel('t, с'); ylabel('h(t)');
    legend(sweeps{j, 3}, 'Location', 'best');
end
