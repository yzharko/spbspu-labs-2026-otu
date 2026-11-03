clear; clc; close all;

Ks = sym(11)/10;
Ts = sym(1)/5;
T1s = sym(3)/10;
T2s = sym(1)/2;
T1as = sym(3)/2;

K = double(Ks);
T = double(Ts);
T1 = double(T1s);
T2 = double(T2s);
T1a = double(T1as);

syms s t

names = {'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
    'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', ...
    'Колебательное', 'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее 1-го порядка', 'Изодромное'};

W = {Ks/s, Ks*s, Ks, Ks/(Ts*s + 1), ...
    Ks/(T2s^2*s^2 + T1as*s + 1), Ks/(T2s^2*s^2 + T1s*s + 1), ...
    Ks/(Ts*s^2 + 1), Ks/(s*(Ts*s + 1)), Ks*s/(Ts*s + 1), ...
    Ks*(Ts*s + 1), Ks*(Ts*s + 1)/s};

Wtf = {tf(K, [1 0]), [], tf(K, 1), tf(K, [T 1]), ...
    tf(K, [T2^2 T1a 1]), tf(K, [T2^2 T1 1]), tf(K, [T 0 1]), ...
    tf(K, [T 1 0]), tf([K 0], [T 1]), [], tf(K*[T 1], [1 0])};

tend = [5 1 1 1.5 8 10 10 5 1.5 1 5];

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

h = cell(1, 11);
w = cell(1, 11);
err = nan(11, 2);

for i = 1:11
    Y_step = W{i}*U_step;
    Y_impulse = W{i}*U_impulse;
    h{i} = ilaplace(Y_step, s, t);
    w{i} = ilaplace(Y_impulse, s, t);

    fprintf('%2d. %s\n', i, names{i});
    fprintf('    h(t) = %s\n', char(h{i}));
    fprintf('    w(t) = %s\n', char(w{i}));

    time = linspace(0, tend(i), 2001);
    [hr, h0, h1] = split_impulses(h{i}, t);
    [wr, w0, w1] = split_impulses(w{i}, t);
    yh = eval_on_grid(hr, t, time);
    yw = eval_on_grid(wr, t, time);

    figure('Name', names{i}, 'Color', 'w');

    subplot(1, 2, 1);
    hold on; grid on;
    plot(time, yh, 'LineWidth', 1.5, 'DisplayName', 'Аналитическое решение');
    draw_impulses(h0, h1, tend(i));
    if ~isempty(Wtf{i})
        ys = step(Wtf{i}, time);
        plot(time(1:100:end), ys(1:100:end), 'ko', 'MarkerSize', 4, 'DisplayName', 'step');
        err(i, 1) = max(abs(ys(:).' - yh));
    end
    title('Переходная функция h(t)');
    xlabel('Время t, с');
    ylabel('h(t)');
    legend('show', 'Location', 'best');
    hold off;

    subplot(1, 2, 2);
    hold on; grid on;
    plot(time, yw, 'LineWidth', 1.5, 'DisplayName', 'Аналитическое решение');
    draw_impulses(w0, w1, tend(i));
    if ~isempty(Wtf{i}) && i ~= 3
        yi = impulse(Wtf{i}, time);
        plot(time(1:100:end), yi(1:100:end), 'ko', 'MarkerSize', 4, 'DisplayName', 'impulse');
        err(i, 2) = max(abs(yi(:).' - yw));
    end
    title('Функция веса w(t)');
    xlabel('Время t, с');
    ylabel('w(t)');
    legend('show', 'Location', 'best');
    hold off;

    sgtitle(names{i});
end

fprintf('\n%-4s %-30s %-14s %-14s\n', 'N', 'Звено', 'max|h - step|', 'max|w - impulse|');
for i = 1:11
    fprintf('%-4d %-30s %-14.2e %-14.2e\n', i, names{i}, err(i, 1), err(i, 2));
end

function [reg, c0, c1] = split_impulses(f, t)
    D0 = sym('D0');
    D1 = sym('D1');
    g = subs(f, [dirac(t), dirac(1, t)], [D0, D1]);
    c0 = double(subs(diff(g, D0), [D0, D1], [0, 0]));
    c1 = double(subs(diff(g, D1), [D0, D1], [0, 0]));
    reg = subs(g, [D0, D1], [0, 0]);
end

function y = eval_on_grid(f, t, time)
    y = real(double(subs(f, t, time))) + zeros(size(time));
end

function draw_impulses(c0, c1, span)
    red = [0.64 0.08 0.18];
    if c0 ~= 0
        stem(0, c0, '^', 'filled', 'Color', red, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%.4g*delta(t)', c0));
    end
    if c1 ~= 0
        dx = 0.035*span;
        stem([-dx dx], [1 -1], 'filled', 'Color', red, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('%.4g*delta''(t)', c1));
    end
end
