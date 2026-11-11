clear; clc; close all;

s = tf("s");

K = 1.7;
T = 0.3;
T1 = 0.4;
T2 = 0.6;
T2a = 0.1;

Ws = {
    withName("Интегрирующее",               @(K, T) K / s)
    withName("Дифференцирующее",            @(K, T) K * s)
    withName("Усилительное",                @(K, T) K)
    withName("Апериодическое 1-го порядка", @(K, T) K / (T * s + 1))
    withName("Апериодическое 2-го порядка", @(K, T) K / ...
                                             (T2a^2 * s^2 + T1 * s + 1))
    withName("Колебательное",               @(K, T) K / ...
                                             (T2^2 * s^2 + T1 * s + 1))
    withName("Консервативное",              @(K, T) K / (T * s^2 + 1))
    withName("Инерционное интегрирующее",   @(K, T) K / (s * (T * s + 1)))
    withName("Реальное дифференцирующее",   @(K, T) K * s / (T * s + 1))
    withName("Форсирующее 1-го порядка",    @(K, T) K * (T * s + 1))
    withName("Изодромное ",                 @(K, T) K * (T * s + 1) / s)
};

coeffs =  [1/5, 1/2, 1, 2, 5];
% Выставляем цвета, чтобы было понятно, в каком направлении идём
% Больше = краснее
t = linspace(0, 1, numel(coeffs))';
red  = [1 0 0];
blue = [0 0 1];
colors = (1 - t) * blue + t * red;
set(groot, 'defaultAxesColorOrder', colors)

% Сначала проварьируем K
batch("K ЛАЧХ", @bodemag,   Ws, K * coeffs, T, 4, 2, 2);
batch("K ЛФЧХ", @bodephase, Ws, K * coeffs, T, 4, 2, 2);
batch("K АФЧХ", @nyquist,   Ws, K * coeffs, T, 4, 2, 2);

% Теперь T
batch("T ЛАЧХ", @bodemag,   Ws, K, T * coeffs, 4, 2, 2);
batch("T ЛФЧХ", @bodephase, Ws, K, T * coeffs, 4, 2, 2);
batch("T АФЧХ", @nyquist,   Ws, K, T * coeffs, 4, 2, 2);

% Бьёт функции на пакеты, чтобы в отчёт всё красиво поместилось
function batch(name, plotFn, Ws, K_range, T_range, bsize, rows, cols)
    for start = 1:bsize:numel(Ws)
        finish = min(start + bsize - 1, numel(Ws));
        subName = strcat(name, " ", string(start), "-", string(finish));
        plotAll(subName, plotFn, Ws(start:finish), ...
            K_range, T_range, rows, cols);
    end
end

% Рисует все переданные функции
function plotAll(name, plotFn, Ws, K_range, T_range, rows, cols)
    figure("Name", name);
    tiledlayout(rows, cols, ...
        'TileSpacing','compact', 'Padding','compact');
    for i = 1:numel(Ws)
        nexttile;
        W = Ws{i};
        sys = cell(numel(K_range) * numel(T_range), 1);
        idx = 0;
        for k = 1:numel(K_range)
            for t = 1:numel(T_range)
                idx = idx + 1;
                sys{idx} = tf(W.expr(K_range(k), T_range(t)));
            end
        end
        plotFn(sys{:});
        grid on;
        title(W.name);
    end
end

function bodephase(varargin)
    bp = bodeplot(varargin{:});
    bp.MagnitudeVisible = "off";
end

function res = withName(name, expr)
    res = struct("name", name, "expr", expr);
end
