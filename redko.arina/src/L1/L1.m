clear; clc; close all;

syms s t K;
T = 0.9; 
T1 = 1.0; 
T2 = 1.4;
T1a = 3.4; 
T2a = 1.4;
s_tf = tf('s');

links = {
    struct('name', 'Интегрирующее', 'W', K/s, 'W_check', @(k) k/s_tf)
    struct('name', 'Идеальное дифференцирующее', 'W', K*s, 'W_check', @(k) k*s_tf)
    struct('name', 'Усилительное', 'W', K, 'W_check', @(k) tf(k))
    struct('name', 'Апериодическое 1-го порядка', 'W', K/(T*s+1), 'W_check', @(k) k/(T*s_tf+1))
    struct('name', 'Апериодическое 2-го порядка', 'W', K/(T2a^2*s^2+T1a*s+1), 'W_check', @(k) k/(T2a^2*s_tf^2+T1a*s_tf+1))
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

K_values = 3.2;
time = 0:0.01:10;

figure;
tiledlayout(4, 3);
for i = 1:N
    nexttile;
    hold on;
    grid on;
    h_k = subs(links{i}.h, K, K_values);
    h_k = subs(h_k, {heaviside(t), dirac(t), dirac(1,t)}, {1, 0, 0});
    y = double(subs(h_k, t, time)) .* ones(size(time));
    plot(time, y);
    title(sprintf('%d. %s', i, links{i}.name));
    xlabel('Время t, с');
    ylabel('h(t)');
    hold off;
end
sgtitle('Переходные характеристики динамических звеньев');

figure;
tiledlayout(4, 3);
for i = 1:N
    nexttile;
    hold on;
    grid on;
    w_k = subs(links{i}.w, K, K_values);
    w_k = subs(w_k, {heaviside(t), dirac(t), dirac(1,t)}, {1, 0, 0});
    y = double(subs(w_k, t, time)) .* ones(size(time));
    plot(time, y);
    title(sprintf('%d. %s', i, links{i}.name));
    xlabel('Время t, с');
    ylabel('w(t)');
    hold off;
end
sgtitle('Импульсные характеристики динамических звеньев');

dt = time(2) - time(1);

figure;
tiledlayout(4, 3);
for i = 1:N
    nexttile;
    hold on;
    grid on;
    if i == 2
        y = zeros(size(time));  y(1) = K_values/dt;
        plot(time, y);
    elseif i == 10
        y = K_values*ones(size(time));  y(1) = K_values + K_values*T/dt;
        plot(time, y);
    else
        [y_step, t_step] = step(links{i}.W_check(K_values), time);
        plot(t_step, y_step, 'DisplayName', sprintf('K=%g', K_values));
    end
    title(sprintf('%d. %s', i, links{i}.name));
    xlabel('Время t, с');
    ylabel('y');
    hold off;
end
sgtitle('Проверка переходной характеристики');

figure;
tiledlayout(4, 3);
for i = 1:N
    nexttile;
    hold on;
    grid on;
    if ismember(i, improper)
        title(sprintf('%d. %s (проверка неприменима)', i, links{i}.name));
    else
        [y, t_out] = impulse(links{i}.W_check(K_values), time);
        plot(t_out, y, 'DisplayName', sprintf('K = %g', K_values));
        title(sprintf('%d. %s', i, links{i}.name));
    end
    xlabel('Время t, с');
    ylabel('y');
    hold off;
end
sgtitle('Проверка импульсной характеристики');

K_var = [0.5 1 3.2 6];
td = 0:0.01:10;  te = 0:0.01:30;
f_aperiodic_K = figure; hold on; grid on;
f_oscillatory_K = figure; hold on; grid on;
for k = K_var
    figure(f_aperiodic_K); plot(td, step(tf(k, [T 1]), td), 'DisplayName', sprintf('K = %g', k));
    figure(f_oscillatory_K); plot(te, step(tf(k, [T2^2 T1 1]), te), 'DisplayName', sprintf('K = %g', k));
end
figure(f_aperiodic_K); title('Звено 4: h(t) при изменении K'); xlabel('Время t, с'); ylabel('h(t)'); legend('show');
figure(f_oscillatory_K); title('Звено 6: h(t) при изменении K'); xlabel('Время t, с'); ylabel('h(t)'); legend('show');

K0 = 3.2;

ta = 0:0.01:25;
f_aperiodic_h_T = figure; hold on; grid on;
f_aperiodic_w_T = figure; hold on; grid on;
fprintf('Звено 4: влияние T при K = %g\n', K0);
for Tk = T*[1/5 1/2 1 2 5]
    W4 = tf(K0, [Tk 1]);
    figure(f_aperiodic_h_T); plot(ta, step(W4, ta),    'DisplayName', sprintf('T = %g', Tk));
    figure(f_aperiodic_w_T); plot(ta, impulse(W4, ta), 'DisplayName', sprintf('T = %g', Tk));
    info = stepinfo(W4, 'SettlingTimeThreshold', 0.05);
    fprintf('T = %-5g tп(5%%) = %5.2f c   tп/T = %.2f   w(0) = K/T = %.2f\n', Tk, info.SettlingTime, info.SettlingTime/Tk, K0/Tk);
end
figure(f_aperiodic_h_T); title('Звено 4: h(t) при изменении T'); xlabel('Время t, с'); ylabel('h(t)'); legend('show');
figure(f_aperiodic_w_T); title('Звено 4: w(t) при изменении T'); xlabel('Время t, с'); ylabel('w(t)'); legend('show');

tc = 0:0.05:300;
f_oscillatory_h_T1 = figure; hold on; grid on;
f_oscillatory_w_T1 = figure; hold on; grid on;
fprintf('Звено 6: влияние T1 при T2 = %g, K = %g\n', T2, K0);
for T1k = [1 0.5 0.25 0.1]
    W6 = tf(K0, [T2^2 T1k 1]);
    figure(f_oscillatory_h_T1); plot(tc, step(W6, tc),    'DisplayName', sprintf('T1 = %g', T1k));
    figure(f_oscillatory_w_T1); plot(tc, impulse(W6, tc), 'DisplayName', sprintf('T1 = %g', T1k));
    info = stepinfo(W6, 'SettlingTimeThreshold', 0.05);
    fprintf('T1 = %-5g zeta = %.3f   перерегулирование = %5.2f %%   tп(5%%) = %6.2f c\n', T1k, T1k/(2*T2), info.Overshoot, info.SettlingTime);
end
figure(f_oscillatory_h_T1); title('Звено 6: h(t) при уменьшении T1'); xlabel('Время t, с'); ylabel('h(t)'); xlim([0 80]); legend('show');
figure(f_oscillatory_w_T1); title('Звено 6: w(t) при уменьшении T1'); xlabel('Время t, с'); ylabel('w(t)'); xlim([0 80]); legend('show');
