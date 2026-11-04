clear; clc; close all;

K0 = 3;
T0 = 1;
T1 = 1.2;
T2 = 1.5;
T1a = 3.75;

syms s t
syms K T positive

names = {'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
    'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', 'Колебательное', ...
    'Консервативное', 'Инерционное интегрирующее', 'Реальное дифференцирующее', ...
    'Форсирующее 1-го порядка', 'Изодромное'};
W = [K/s; K*s; K; K/(T*s + 1); K/(T2^2*s^2 + T1a*s + 1); K/(T2^2*s^2 + T1*s + 1); ...
    K/(T*s^2 + 1); K/(s*(T*s + 1)); K*s/(T*s + 1); K*(T*s + 1); K*(T*s + 1)/s];
N = numel(names);

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

h = sym(zeros(N, 1));
w = sym(zeros(N, 1));
for i = 1:N
    Y_step = W(i)*U_step;
    Y_impulse = W(i)*U_impulse;
    h(i) = simplify(ilaplace(Y_step, s, t));
    w(i) = simplify(ilaplace(Y_impulse, s, t));
    fprintf('%d. %s\n   h(t) = %s\n   w(t) = %s\n', i, names{i}, char(h(i)), char(w(i)));
end

R = subs([h w], [K T], [K0 T0]);
R0 = subs(subs(R, dirac(1, t), 0), dirac(t), 0);
D0 = double(subs(subs(R, dirac(1, t), 0), dirac(t), 1) - R0);
D1 = double(subs(subs(R, dirac(t), 0), dirac(1, t), 1) - R0);

s_tf = tf('s');
W_check = {K0/s_tf, [], tf(K0), K0/(T0*s_tf + 1), K0/(T2^2*s_tf^2 + T1a*s_tf + 1), ...
    K0/(T2^2*s_tf^2 + T1*s_tf + 1), K0/(T0*s_tf^2 + 1), K0/(s_tf*(T0*s_tf + 1)), ...
    K0*s_tf/(T0*s_tf + 1), [], K0*(T0*s_tf + 1)/s_tf};

time = 0:0.01:15;
simulate = {@step, @impulse};
y_labels = {'h(t)', 'w(t)'};
fig_titles = {'Переходные характеристики h(t)', 'Импульсные характеристики w(t)'};
dev = NaN(N, 2);
for k = 1:2
    figure('Name', fig_titles{k}, 'Position', [50 50 1000 1080]);
    tiledlayout(ceil(N/3), 3);
    for i = 1:N
        nexttile; hold on; grid on;
        y = double(subs(R0(i, k), t, time)) + zeros(size(time));
        plot([-1 0 time], [0 0 y], 'LineWidth', 3);
        if ~isempty(W_check{i})
            y_sim = simulate{k}(W_check{i}, time);
            plot(time, y_sim, '--', 'LineWidth', 1.2);
            dev(i, k) = max(abs(y_sim' - y));
        end
        if i == 1
            legend('аналитический расчёт', 'моделирование', 'Location', 'best');
        end
        if D0(i, k) ~= 0
            stem(0, D0(i, k), 'r^', 'filled');
        end
        if D1(i, k) ~= 0
            stem(0.5, D1(i, k), 'r^', 'filled');
            stem(0.5, -D1(i, k), 'rv', 'filled');
        end
        title(sprintf('%d. %s', i, names{i}));
        xlabel('Время t, с');
        ylabel(y_labels{k});
        xlim([-1 time(end)]);
    end
    sgtitle(sprintf('%s при K = %g, T = %g с', fig_titles{k}, K0, T0));
end
fprintf('\nНаибольшее расхождение аналитики и моделирования: %.1e\n', max(dev(:)));

figure('Name', 'Влияние K', 'Position', [50 50 1000 400]);
tiledlayout(1, 2);
ax4 = nexttile; hold on; grid on;
ax6 = nexttile; hold on; grid on;
for Ki = K0*[0.5 1 2]
    plot(ax4, time, step(Ki/(T0*s_tf + 1), time), 'DisplayName', sprintf('K = %g', Ki));
    plot(ax6, time, step(Ki/(T2^2*s_tf^2 + T1*s_tf + 1), time), ...
        'DisplayName', sprintf('K = %g', Ki));
end
title(ax4, '4. Апериодическое 1-го порядка');
title(ax6, '6. Колебательное');
for ax = [ax4 ax6]
    xlabel(ax, 'Время t, с');
    ylabel(ax, 'h(t)');
    legend(ax, 'Location', 'southeast');
end

time_aperiodic = 0:0.01:30;
figure('Name', 'Апериодическое звено', 'Position', [50 50 1000 400]);
tiledlayout(1, 2);
ax_h = nexttile; hold on; grid on;
ax_w = nexttile; hold on; grid on;
fprintf('\nАпериодическое звено: T, tп(5%%), tн(10-90%%), w(0)\n');
for Ti = T0*[0.2 0.5 1 2 5]
    W4 = K0/(Ti*s_tf + 1);
    y_h = step(W4, time_aperiodic);
    y_w = impulse(W4, time_aperiodic);
    info = stepinfo(y_h, time_aperiodic, K0, 'SettlingTimeThreshold', 0.05, ...
        'RiseTimeLimits', [0.1 0.9]);
    plot(ax_h, time_aperiodic, y_h, 'DisplayName', sprintf('T = %g с', Ti));
    plot(ax_w, time_aperiodic, y_w, 'DisplayName', sprintf('T = %g с', Ti));
    fprintf('%5g %7.2f %6.2f %6.2f\n', Ti, info.SettlingTime, info.RiseTime, y_w(1));
end
sgtitle(sprintf('Апериодическое звено при K = %g и различных T', K0));
for ax = [ax_h ax_w]
    xlabel(ax, 'Время t, с');
    legend(ax, 'Location', 'east');
end
ylabel(ax_h, 'h(t)');
ylabel(ax_w, 'w(t)');

time_oscillatory = 0:0.01:200;
figure('Name', 'Колебательное звено', 'Position', [50 50 1000 700]);
tiledlayout(2, 1);
ax_h = nexttile; hold on; grid on;
ax_w = nexttile; hold on; grid on;
fprintf('\nКолебательное звено: T1, σ %%, tп(5%%), период\n');
for T1i = T1*[1 0.5 0.25 0.1 0]
    W6 = K0/(T2^2*s_tf^2 + T1i*s_tf + 1);
    y_h = step(W6, time_oscillatory);
    y_w = impulse(W6, time_oscillatory);
    info = stepinfo(y_h, time_oscillatory, K0, 'SettlingTimeThreshold', 0.05);
    period = 2*pi*2*T2^2/sqrt(4*T2^2 - T1i^2);
    plot(ax_h, time_oscillatory, y_h, 'DisplayName', sprintf('T1 = %g с', T1i));
    plot(ax_w, time_oscillatory, y_w, 'DisplayName', sprintf('T1 = %g с', T1i));
    fprintf('%5g %7.2f %7.2f %6.2f\n', T1i, info.Overshoot, info.SettlingTime, period);
end
sgtitle(sprintf('Колебательное звено при K = %g, T2 = %g с и уменьшении T1', K0, T2));
for ax = [ax_h ax_w]
    xlim(ax, [0 80]);
    xlabel(ax, 'Время t, с');
    legend(ax, 'Location', 'eastoutside');
end
ylabel(ax_h, 'h(t)');
ylabel(ax_w, 'w(t)');
