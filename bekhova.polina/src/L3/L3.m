clc;
clear;
close all;


K = 3.8;
T = 1.1;
T1 = 1.5;
T2 = 1.6;
T1_ap = 3.5;

s = tf('s');
w = logspace(-2, 2, 1000);

K_values = K * [0.25 0.5 1 2 4];
T_values = T * [0.25 0.5 1 2 4];

T_grad = [0.22, 0.55, 1.1, 2.2, 5.5];

names = {'Интегрирующее', 'Идеальное дифференцирующее', ...
    'Усилительное', 'Апериодическое 1-го порядка', ...
    'Апериодическое 2-го порядка', 'Колебательное', ...
    'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее', 'Изодромное'};

N = numel(names);


W = cell(1, N);
for n = 1:N
    W{n} = make_tf(n, K, T, T1, T2, T1_ap, s);
end

plotNames = {'АФЧХ', 'ЛАЧХ', 'ЛФЧХ'};

for p = 1:3
    figure(p);
    set(gcf, 'Name', plotNames{p}, 'NumberTitle', 'off', ...
        'Position', [50 50 1500 850]);
    tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    for n = 1:N
        ax = nexttile;
        plot_characteristic(ax, W{n}, w, p);
        title(ax, sprintf('%d. %s', n, names{n}), 'FontSize', 8);
    end
    sgtitle([plotNames{p} ' всех динамических звеньев']);
end


for p = 1:2
    figure(3 + p);
    set(gcf, 'Name', 'Влияние коэффициента K', 'NumberTitle', 'off', ...
        'Position', [50 50 1500 850]);
    tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    for n = 1:N
        ax = nexttile;
        hold(ax, 'on');

        for i = 1:numel(K_values)
            Wvar = make_tf(n, K_values(i), T, T1, T2, T1_ap, s);
            plot_characteristic(ax, Wvar, w, p + 1);
        end

        title(ax, sprintf('%d. %s', n, names{n}), 'FontSize', 8);
        if n == N
            legend(ax, compose('K = %.2g', K_values), ...
                'Location', 'best', 'FontSize', 7);
        end
    end
    sgtitle(['Влияние K на ' plotNames{p + 1}]);
end


for p = 1:2
    figure(5 + p);
    set(gcf, 'Name', 'Влияние постоянных времени', ...
        'NumberTitle', 'off', 'Position', [50 50 1500 850]);
    tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    for n = 1:N
        ax = nexttile;
        hold(ax, 'on');

        if n <= 3
            plot_characteristic(ax, W{n}, w, p + 1);
        else
            for i = 1:numel(T_values)
                Tv = T_values(i);

                if n == 5
                    Wvar = K / (Tv^2 * s^2 + T1_ap * Tv / T2 * s + 1);
                elseif n == 6
                    Wvar = K / (Tv^2 * s^2 + T1 * Tv / T2 * s + 1);
                else
                    Wvar = make_tf(n, K, Tv, T1, T2, T1_ap, s);
                end

                plot_characteristic(ax, Wvar, w, p + 1);
            end
        end

        title(ax, sprintf('%d. %s', n, names{n}), 'FontSize', 8);
        if n == N
            legend(ax, compose('T = %.2g', T_values), ...
                'Location', 'best', 'FontSize', 7);
        end
    end
    sgtitle(['Влияние постоянных времени на ' plotNames{p + 1}]);
end


figure(8);
set(gcf, 'Name', 'АФЧХ звена №4', 'NumberTitle', 'off', ...
    'Position', [100 100 1250 550]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile;
hold(ax, 'on');
for i = 1:numel(K_values)
    Wvar = K_values(i) / (T * s + 1);
    plot_characteristic(ax, Wvar, w, 1);
end
title(ax, 'Изменение коэффициента K');
legend(ax, compose('K = %.2g', K_values), 'Location', 'best');

ax = nexttile;
hold(ax, 'on');
for i = 1:numel(T_values)
    Wvar = K / (T_values(i) * s + 1);
    plot_characteristic(ax, Wvar, w, 1);
end
title(ax, 'Изменение постоянной времени T');
legend(ax, compose('T = %.2g', T_values), 'Location', 'best');

sgtitle('Сравнение АФЧХ апериодического звена первого порядка');



figure(9);
set(gcf, 'Name', 'Исследование звена №5', 'NumberTitle', 'off', ...
    'Position', [100 100 1250 550]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

T1_values = [3.2, 3.5, 4.0, 5.0, 6.4];
ax = nexttile;
hold(ax, 'on');

for i = 1:numel(T1_values)
    Wvar = K / (T2^2 * s^2 + T1_values(i) * s + 1);
    plot_characteristic(ax, Wvar, w, 2);
end

title(ax, 'Изменение T_{1ап} при постоянном T_2');
legend(ax, compose('T_{1ап} = %.2g', T1_values), 'Location', 'best');

T2_values_ap = [0.22, 0.55, 1.1, 1.5, 1.75];
ax = nexttile;
hold(ax, 'on');

for i = 1:numel(T2_values_ap)
    Wvar = K / (T2_values_ap(i)^2 * s^2 + T1_ap * s + 1);
    plot_characteristic(ax, Wvar, w, 2);
end

title(ax, 'Изменение T_2 при постоянном T_{1ап}');
legend(ax, compose('T_2 = %.2g', T2_values_ap), 'Location', 'best');

sgtitle('Исследование апериодического звена второго порядка');


figure(10);
set(gcf, 'Name', 'Исследование звена №6', 'NumberTitle', 'off', ...
    'Position', [100 100 1250 550]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

T1_values_osc = T_grad(T_grad < 2*T2);

ax = nexttile;
hold(ax, 'on');
for i = 1:numel(T1_values_osc)
    Wvar = K / (T2^2 * s^2 + T1_values_osc(i) * s + 1);
    plot_characteristic(ax, Wvar, w, 2);
end
title(ax, 'Изменение T_1 при постоянном T_2');
legend(ax, compose('T_1 = %.2g', T1_values_osc), 'Location', 'best');

T2_values_osc = T_grad(T_grad > T1/2);

ax = nexttile;
hold(ax, 'on');
for i = 1:numel(T2_values_osc)
    Wvar = K / (T2_values_osc(i)^2 * s^2 + T1 * s + 1);
    plot_characteristic(ax, Wvar, w, 2);
end
title(ax, 'Изменение T_2 при постоянном T_1');
legend(ax, compose('T_2 = %.2g', T2_values_osc), 'Location', 'best');

sgtitle('Исследование колебательного звена второго порядка');


function W = make_tf(n, K, T, T1, T2, T1_ap, s)
    switch n
        case 1,  W = K / s;
        case 2,  W = K * s;
        case 3,  W = tf(K);
        case 4,  W = K / (T * s + 1);
        case 5,  W = K / (T2^2 * s^2 + T1_ap * s + 1);
        case 6,  W = K / (T2^2 * s^2 + T1 * s + 1);
        case 7,  W = K / (T * s^2 + 1);
        case 8,  W = K / (s * (T * s + 1));
        case 9,  W = K * s / (T * s + 1);
        case 10, W = K * (T * s + 1);
        case 11, W = K * (T * s + 1) / s;
    end
end

function plot_characteristic(ax, W, w, type)
    persistent colors
    if isempty(colors)
        colors = lines(5);
    end

    if type == 1
        H = squeeze(freqresp(W, w));
        H = H(:);
        plot(ax, real(H), imag(H), 'LineWidth', 1.2);
        hold(ax, 'on');
        plot(ax, real(H), -imag(H), 'LineWidth', 1.2);
        xline(ax, 0, 'k:');
        yline(ax, 0, 'k:');
        xlabel(ax, 'Вещественная часть');
        ylabel(ax, 'Мнимая часть');

    elseif type == 2
        [mag, ~] = bode(W, w);
        semilogx(ax, w, 20 * log10(squeeze(mag)), ...
            'LineWidth', 1.2);
        xlabel(ax, 'Частота, рад/с');
        ylabel(ax, 'Амплитуда, дБ');

    else
        [~, phase] = bode(W, w);
        semilogx(ax, w, squeeze(phase), 'LineWidth', 1.2);
        xlabel(ax, 'Частота, рад/с');
        ylabel(ax, 'Фаза, град');
    end

    grid(ax, 'on');
    box(ax, 'on');
    set(ax, 'FontSize', 8);
end
