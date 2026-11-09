clear;
clc;
close all;

K = 1.8;
T = 0.6;
T1 = 0.9;
T2 = 1.2;
T1_aper = 2.6; % Для 5 функции T1_aper >= 2*T2

if T1_aper < 2*T2
    error('Для 5-го звена необходимо T1_aper >= 2*T2.');
end

time = 0:0.01:20;

names = { ...
    'Интегрирующее (идеальное)', ...
    'Дифференцирующее (идеальное)', ...
    'Усилительное (безынерционное)', ...
    'Апериодическое 1-го порядка', ...
    'Апериодическое 2-го порядка', ...
    'Колебательное', ...
    'Консервативное', ...
    'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', ...
    'Форсирующее 1-го порядка', ...
    'Изодромное'};

%% Исследование всех 11 типовых динамических звеньев
syms s t
syms K_s T_s T1_s T2_s T1_aper_s positive 
% positive - ключевое слово матлаба, говорит что эти символьные параметры
% положительные и матлабу легче упрощать имея это ввиду. К параметрам s и t
% это не относится, поэтому декларация разбита на 2 строчки

W_formula = cell(1, 11);
W_formula{1}  = K_s/s;
W_formula{2}  = K_s*s;
W_formula{3}  = K_s;
W_formula{4}  = K_s/(T_s*s + 1);
W_formula{5}  = K_s/(T2_s^2*s^2 + T1_aper_s*s + 1);
W_formula{6}  = K_s/(T2_s^2*s^2 + T1_s*s + 1);
W_formula{7}  = K_s/(T_s^2*s^2 + 1);
W_formula{8}  = K_s/(s*(T_s*s + 1));
W_formula{9}  = K_s*s/(T_s*s + 1);
W_formula{10} = K_s*(T_s*s + 1);
W_formula{11} = K_s*(T_s*s + 1)/s;

parameter_symbols = [K_s, T_s, T1_s, T2_s, T1_aper_s];
parameter_values = [K, T, T1, T2, T1_aper];

U_step = laplace(heaviside(t), t, s);
U_impulse = laplace(dirac(t), t, s);

W_sym = cell(1, 11);
Y_step_sym = cell(1, 11);
h_sym = cell(1, 11);
Y_impulse_sym = cell(1, 11);
w_sym = cell(1, 11);

Y_step_regular_sym = cell(1, 11);
h_regular_sym = cell(1, 11);
Y_impulse_regular_sym = cell(1, 11);
w_regular_sym = cell(1, 11);

for i = 1:numel(names)
    W_sym{i} = simplify(subs(W_formula{i}, parameter_symbols, parameter_values));

    Y_step_sym{i} = simplify(W_sym{i} * U_step);
    h_sym{i} = simplify(ilaplace(Y_step_sym{i}, s, t));

    Y_impulse_sym{i} = simplify(W_sym{i} * U_impulse);
    w_sym{i} = simplify(ilaplace(Y_impulse_sym{i}, s, t));

    [Y_step_regular_sym{i}, h_regular_sym{i}] = ...
        finitePartFromLaplace(Y_step_sym{i}, s, t);

    [Y_impulse_regular_sym{i}, w_regular_sym{i}] = ...
        finitePartFromLaplace(Y_impulse_sym{i}, s, t);
    
    fprintf('Звено %d: %s\n', i, names{i});
    fprintf('W(s) = %s\n', char(W_sym{i}));
    fprintf('h(t) = %s\n', char(h_sym{i}));
    fprintf('w(t) = %s\n', char(w_sym{i}));
    fprintf('Конечная часть h(t) для построения графика = %s\n', char(h_regular_sym{i}));
    fprintf('Конечная часть w(t) для построения графика = %s\n\n', char(w_regular_sym{i}));
end

h_fun = cell(1, 11);
w_fun = cell(1, 11);
for i = 1:11
    h_fun{i} = matlabFunction(h_regular_sym{i}, 'Vars', t);
    w_fun{i} = matlabFunction(w_regular_sym{i}, 'Vars', t);
end

%% Построение графиков h(t) и w(t)
plotResponses(time, h_fun, names, ...
    'Ступенчатое входное воздействие', 'h(t)', 'Время (s)', 'h(t)', 'subplots', 1.3);

plotResponses(time, w_fun, names, ...
    'Импульсное входное воздействие', 'w(t)', 'Время (s)', 'w(t)', 'subplots', 1.3);

%% Определение влияния варьирования коэффициента K на характер переходного процесса на примере 6-го звена
K_values = 1:10;
h_K_fun = cell(size(K_values));
w_K_fun = cell(size(K_values));

for j = 1:numel(K_values)
    values_now = [K_values(j), T, T1, T2, T1_aper];
    W_variant = simplify(subs(W_formula{6}, parameter_symbols, values_now));

    [~, h_variant] = finitePartFromLaplace(W_variant/s, s, t);
    h_K_fun{j} = matlabFunction(h_variant, 'Vars', t);

    [~, w_variant] = finitePartFromLaplace(W_variant, s, t);
    w_K_fun{j} = matlabFunction(w_variant, 'Vars', t);
end

K_labels = arrayfun(@(value) sprintf('K = %g', value), ...
    K_values, 'UniformOutput', false);

plotResponses(time, h_K_fun, K_labels, ...
    'Влияние изменения K при ступенчатом воздействии на 6 звено', ...
    'Влияние изменения K при ступенчатом воздействии на 6 звено', ...
    'Время (s)', 'h(t)', 'overlay', 1.2);

plotResponses(time, w_K_fun, K_labels, ...
    'Влияние изменения K при импульсном воздействии на 6 звено', ...
    'Влияние изменения K при импульсном воздействии на 6 звено', ...
    'Время (s)', 'w(t)', 'overlay', 1.2);

%% Изменение скорости протекания процесса для 4 звена при увеличении и уменьшении постоянной времени T в 2 и в 5 раз
T_values = T * [1/5 1/2 1 2 5];
h_T_fun = cell(size(T_values));
w_T_fun = cell(size(T_values));

for j = 1:numel(T_values)
    values_now = [K, T_values(j), T1, T2, T1_aper];
    W_variant = simplify(subs(W_formula{4}, parameter_symbols, values_now));

    [~, h_variant] = finitePartFromLaplace(W_variant/s, s, t);
    h_T_fun{j} = matlabFunction(h_variant, 'Vars', t);

    [~, w_variant] = finitePartFromLaplace(W_variant, s, t);
    w_T_fun{j} = matlabFunction(w_variant, 'Vars', t);
end

T_labels = arrayfun(@(value) sprintf('T = %.3g s', value), ...
    T_values, 'UniformOutput', false);

plotResponses(time, h_T_fun, T_labels, ...
    'Влияние изменения T при ступенчатом воздействии на 4 звено', ...
    'Влияние изменения T при ступенчатом воздействии на 4 звено', ...
    'Время (s)', 'h(t)', 'overlay', 1.2);

plotResponses(time, w_T_fun, T_labels, ...
    'Влияние изменения T при импульсном воздействии на 4 звено', ...
    'Влияние изменения T при импульсном воздействии на 4 звено', ...
    'Время (s)', 'w(t)', 'overlay', 1.2);

%% Исследование изменения перерегулирования и времени затухания при уменьшении коэффициента T1 на примере 6-го звена
T1_values = [T1, 2*T1/3, T1/3, T1/10, T1/20];
time_long = 0:0.02:250;
h_T1_fun = cell(size(T1_values));
w_T1_fun = cell(size(T1_values));

for j = 1:numel(T1_values)
    values_now = [K, T, T1_values(j), T2, T1_aper];
    W_variant = simplify(subs(W_formula{6}, parameter_symbols, values_now));

    [~, h_variant] = finitePartFromLaplace(W_variant/s, s, t);
    h_T1_fun{j} = matlabFunction(h_variant, 'Vars', t);

    [~, w_variant] = finitePartFromLaplace(W_variant, s, t);
    w_T1_fun{j} = matlabFunction(w_variant, 'Vars', t);
end

T1_labels = arrayfun(@(value) sprintf('T1 = %.3g s', value), ...
    T1_values, 'UniformOutput', false);

plotResponses(time_long, h_T1_fun, T1_labels, ...
    'Влияние изменения T1 при ступенчатом воздействии на 6 звено', ...
    'Влияние изменения T1 при ступенчатом воздействии на 6 звено', ...
    'Время (s)', 'h(t)', 'overlay', 1.1);

plotResponses(time_long, w_T1_fun, T1_labels, ...
    'Влияние изменения T1 при импульсном воздействии на 6 звено', ...
    'Влияние изменения T1 при импульсном воздействии на 6 звено', ...
    'Время (s)', 'w(t)', 'overlay', 1.1);

fprintf('Изменение перерегулирования и времени затухания при уменьшении коэффициента T1\n');
fprintf(' T1 (s)     zeta       Перерегулирование (%%)    Время затухания (s)\n');
for T1v = T1_values
    zeta = T1v/(2*T2);
    alpha = T1v/(2*T2^2);
    Mp = 100*exp(-zeta*pi/sqrt(1-zeta^2));
    Ts_2pct = 4/alpha;
    fprintf('%8.4f   %8.4f      %8.2f          %10.2f\n', ...
        T1v, zeta, Mp, Ts_2pct);
end

function plotResponses(time, responseFunctions, labels, figureName, ...
    plotTitle, xLabel, yLabel, plotMode, lineWidth)

    figure('Name', figureName, 'Color', 'w');

    switch plotMode
        case 'subplots'
            tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

            for k = 1:numel(responseFunctions)
                nexttile;
                y = responseFunctions{k}(time);

                if isscalar(y)
                    y = y .* ones(size(time));
                end
                y = reshape(y, size(time));

                plot(time, y, 'LineWidth', lineWidth);
                grid on;
                title(sprintf('%d. %s', k, labels{k}), ...
                    'Interpreter', 'none');
                xlabel(xLabel);
                ylabel(yLabel);
            end
            sgtitle(plotTitle);

        case 'overlay'
            hold on;
            grid on;

            for k = 1:numel(responseFunctions)
                y = responseFunctions{k}(time);

                if isscalar(y)
                    y = y .* ones(size(time));
                end
                y = reshape(y, size(time));

                plot(time, y, 'LineWidth', lineWidth, ...
                    'DisplayName', labels{k});
            end

            xlabel(xLabel);
            ylabel(yLabel);
            title(plotTitle, 'Interpreter', 'none');
            legend('show', 'Location', 'best');
            hold off;

        otherwise
            error('Неизвестный режим построения графика: %s', plotMode);
    end
end

function [Y_regular, response_regular] = finitePartFromLaplace(Y, s, t)
    [numerator, denominator] = numden(Y);
    [~, remainder] = quorem(expand(numerator), expand(denominator), s);

    Y_regular = simplify(remainder/denominator);
    response_regular = simplify(ilaplace(Y_regular, s, t));
    response_regular = simplify(subs(response_regular, heaviside(t), sym(1)));
end
