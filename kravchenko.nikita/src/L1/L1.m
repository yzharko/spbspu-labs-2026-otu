clear;
clc;
close all;

K = 1.2;
T = 0.3;
T1 = 0.5;
T2 = 0.6;
T1_5 = 3*T2;

names = {
    'Идеальное интегрирующее'
    'Идеальное дифференцирующее'
    'Усилительное'
    'Апериодическое 1-го порядка'
    'Апериодическое 2-го порядка'
    'Колебательное'
    'Консервативное'
    'Инерционное интегрирующее'
    'Реальное дифференцирующее'
    'Форсирующее 1-го порядка'
    'Изодромное'
};

syms s t real
Ustep = laplace(heaviside(t), t, s);
Uimpulse = laplace(dirac(t), t, s);
fprintf('Вариант  K = %g, T = %g, T1 = %g, T2 = %g\n', K, T, T1, T2);
fprintf('Для звена №5  T1 = %g с > 2T2 = %g с.\n', T1_5, 2*T2);
fprintf('U скачка(s) = %s, U импульса(s) = %s\n', char(Ustep), char(Uimpulse));
fprintf('Отклики причинные. Формулы регулярных частей приведены для t ≥ 0.\n');

[numerators, denominators] = linkPolynomials(K, T, T1, T2, T1_5);
time = (0:0.01:12).';
plotTime = [-0.4; 0; time];
gains = K*[0.5 2 1];
colors = [0.18 0.54 0.38; 0.83 0.37 0.12; 0.10 0.30 0.68];
styles = {'--', '-.', '-'};
responseLabels = {'h(t)', 'w(t)'};
windowTitles = {'Переходные характеристики', 'Импульсные характеристики'};
layouts = gobjects(1, 2);
for response = 1:2
    layouts(response) = makeLayout(windowTitles{response}, 6, 2);
end

for n = 1:numel(names)
    num = numerators{n};
    den = denominators{n};
    W = poly2sym(sym(num, 'r'), s)/poly2sym(sym(den, 'r'), s);
    Y = simplify([W*Ustep, W*Uimpulse]);
    exact = [ilaplace(Y(1), s, t), ilaplace(Y(2), s, t)];
    [numeric, delta] = numericResponses(num, den, time);
    singular = sym(delta, 'r')*[dirac(t); dirac(1, t)];

    fprintf('\n№%d. %s\n', n, names{n});
    fprintf('W(s) = %s\n', char(W));
    fprintf('Y скачка(s) = %s\n', char(Y(1)));
    fprintf('Y импульса(s) = %s\n', char(Y(2)));
    fprintf('h(t) = %s\nw(t) = %s\n', char(exact(1)), char(exact(2)));
    errors = zeros(1, 2);

    for response = 1:2
        regular = subs(exact(response), [dirac(t), dirac(1, t)], [0, 0]);
        assert(isequal(simplify(exact(response) - regular - singular(response)), sym(0)), ...
            'Не совпали импульсные слагаемые звена №%d.', n);
        evaluate = matlabFunction(subs(regular, heaviside(t), 1), 'Vars', t);
        values = evaluate(time) + zeros(size(time));
        errors(response) = max(abs(values - numeric(:, response)));
        assert(errors(response) < 1e-7, 'Не совпали регулярные части звена №%d.', n);

        ax = nexttile(layouts(response), n);
        hold(ax, 'on');
        for j = 1:numel(gains)
            plot(ax, plotTime, [0; 0; values*gains(j)/K], styles{j}, ...
                'Color', colors(j,:), 'LineWidth', 1.4);
        end
        formatAxes(ax, responseLabels{response});
        caption = sprintf('№%d. %s', n, names{n});
        if any(delta(response,:))
            caption = {caption, [deltaLabel(delta(response,:)) ' при базовом K']};
        end
        title(ax, caption, 'FontSize', 9);
    end

    fprintf('Максимальные ошибки сверки со step и impulse  h = %.2e, w = %.2e\n', errors);
end

for response = 1:2
    ax = nexttile(layouts(response), 12);
    hold(ax, 'on');
    for j = 1:numel(gains)
        plot(ax, NaN, NaN, styles{j}, 'Color', colors(j,:), 'LineWidth', 1.4, ...
            'DisplayName', sprintf('K = %g', gains(j)));
    end
    axis(ax, 'off');
    legend(ax, 'Location', 'north', 'Box', 'off');
    text(ax, 0.5, 0.35, {'Линии показывают регулярную часть.', ...
        'δ(t) и δ′(t) при t = 0 указаны над графиком.', ...
        'При изменении K их коэффициенты', 'меняются пропорционально K.'}, ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 9);
end

band = 0.05;
Tvalues = T*[1/5 1/2 1 2 5];
time4 = linspace(0, 5*max(Tvalues), 1501).';
T1values = T1*[1 1/2 1/5 1/10];
alpha = T1values/(2*T2^2);
zeta = T1values/(2*T2);
omegaD = sqrt(1 - zeta.^2)/T2;
envelopeTime = log(1./(band*sqrt(1 - zeta.^2)))./alpha;
time6 = (0:0.01:ceil(max(envelopeTime) + 1)).';

experiment = makeLayout('Звенья №4 и №6 — скорость и затухание', 4, 1);
ax4h = nexttile(experiment, 1);
ax4w = nexttile(experiment, 2);
ax6h = nexttile(experiment, 3);
ax6w = nexttile(experiment, 4);
for ax = [ax4h ax4w ax6h ax6w]
    hold(ax, 'on');
end

fprintf('\n№4. Изменение T при K = %g\n', K);
fprintf('%8s %14s %15s %14s\n', 'T, с', 't63,2%, с', 't5%, с', 'h''(0+)=w(0+)');
for tau = Tvalues
    response4 = numericResponses(K, [tau 1], time4);
    label = sprintf('T = %g с', tau);
    plot(ax4h, time4, response4(:,1), 'LineWidth', 1.3, 'DisplayName', label);
    plot(ax4w, time4, response4(:,2), 'LineWidth', 1.3, 'DisplayName', label);
    fprintf('%8.3f %14.3f %15.3f %14.3f\n', tau, tau, -tau*log(band), K/tau);
end
yline(ax4h, (1-band)*K, ':', '95% от K', 'HandleVisibility', 'off');

fprintf('\n№6. Уменьшение T1 при T2 = %g с и K = %g\n', T2, K);
fprintf('t5%% — момент последнего входа в полосу ±5%% от установившегося уровня.\n');
fprintf('%8s %10s %13s %14s %13s %16s\n', 'T1, с', 'ζ', 'σ, %', 't пика, с', 't5%, с', 't огибающей, с');
for j = 1:numel(T1values)
    a = alpha(j);
    wd = omegaD(j);
    relativeError = @(x) exp(-a*x).*(cos(wd*x) + (a/wd)*sin(wd*x));
    lastPeak = floor(-log(band)*wd/(a*pi));
    peakTime = lastPeak*pi/wd;
    nextZero = peakTime + (pi/2 + atan(a/wd))/wd;
    settlingTime = fzero(@(x) abs(relativeError(x)) - band, [peakTime nextZero]);
    overshoot = 100*exp(-a*pi/wd);
    response6 = numericResponses(K, [T2^2 T1values(j) 1], time6);
    label = sprintf('T1 = %g с', T1values(j));
    plot(ax6h, time6, response6(:,1), 'LineWidth', 1.3, 'DisplayName', label);
    plot(ax6w, time6, response6(:,2), 'LineWidth', 1.3, 'DisplayName', label);
    fprintf('%8.3f %10.4f %13.2f %14.3f %13.3f %16.3f\n', ...
        T1values(j), zeta(j), overshoot, pi/wd, settlingTime, envelopeTime(j));
end
limitResponse = numericResponses(K, [T2^2 0 1], time6);
plot(ax6h, time6, limitResponse(:,1), ':k', 'DisplayName', 'T1 = 0 (предел)');
plot(ax6w, time6, limitResponse(:,2), ':k', 'DisplayName', 'T1 = 0 (предел)');
yline(ax6h, K*[1-band 1+band], '--', 'HandleVisibility', 'off');

experimentAxes = [ax4h ax4w ax6h ax6w];
experimentTitles = {'№4. Переходная функция', '№4. Функция веса', ...
    '№6. Переходная функция', '№6. Функция веса'};
for j = 1:numel(experimentAxes)
    ax = experimentAxes(j);
    formatAxes(ax, responseLabels{1 + mod(j-1, 2)});
    title(ax, experimentTitles{j});
    legend(ax, 'Location', 'best', 'Box', 'off');
end

baseline = [T T1 T2 T1_5];
sweepLinks = [5 5 6 7 8 9 10 11];
sweepColumns = [4 3 3 1 1 1 1 1];
sweepNames = {'T', 'T1', 'T2', 'T1'};
factors = [0.8 1 1.2; 0.8 1 1.2; 0.5 1 1.5; repmat([0.5 1 2], 5, 1)];
sweepLayouts = [makeLayout('Влияние постоянных времени на h(t)', 4, 2), ...
    makeLayout('Влияние постоянных времени на w(t)', 4, 2)];

for j = 1:numel(sweepLinks)
    n = sweepLinks(j);
    column = sweepColumns(j);
    parameter = sweepNames{column};
    axesPair = [nexttile(sweepLayouts(1), j), nexttile(sweepLayouts(2), j)];
    hold(axesPair(1), 'on');
    hold(axesPair(2), 'on');

    for factor = factors(j,:)
        varied = baseline;
        varied(column) = baseline(column)*factor;
        [num, den] = linkPolynomials(K, varied(1), varied(2), varied(3), varied(4));
        [curves, delta] = numericResponses(num{n}, den{n}, time);

        for response = 1:2
            label = sprintf('%s = %g', parameter, varied(column));
            if any(delta(response,:))
                label = [label '  |  ' deltaLabel(delta(response,:))];
            end
            plot(axesPair(response), plotTime, [0; 0; curves(:,response)], ...
                'LineWidth', 1.3, 'DisplayName', label);
        end
    end

    for response = 1:2
        ax = axesPair(response);
        formatAxes(ax, responseLabels{response});
        title(ax, {sprintf('№%d. %s', n, names{n}), ['Меняется ' parameter]}, 'FontSize', 9);
        legend(ax, 'Location', 'best', 'FontSize', 8, 'Box', 'off');
    end
end

function [num, den] = linkPolynomials(K, T, T1, T2, T1_5)
    num = {K; [K 0]; K; K; K; K; K; K; [K 0]; [K*T K]; [K*T K]};
    den = {[1 0]; 1; 1; [T 1]; [T2^2 T1_5 1]; [T2^2 T1 1]; ...
        [T 0 1]; [T 1 0]; [T 1]; 1; [1 0]};
end

function [curves, delta] = numericResponses(num, den, time)
    direct = 0;
    derivative = 0;
    remainder = num;
    if numel(num) >= numel(den)
        [quotient, remainder] = deconv(num, den);
        direct = quotient(end);
        if numel(quotient) == 2
            derivative = quotient(1);
        end
        remainder = remainder(end-numel(den)+2:end);
    end

    if isempty(remainder) || all(remainder == 0)
        curves = [direct*ones(size(time)), zeros(size(time))];
    else
        dynamics = tf(remainder, den);
        curves = [step(dynamics, time) + direct, impulse(dynamics, time)];
    end
    delta = [derivative 0; direct derivative];
end

function layout = makeLayout(caption, rows, columns)
    screenSize = get(groot, 'ScreenSize');
    screenRatio = screenSize(3)/screenSize(4);
    portraitRatio = 0.78;
    height = min(0.83, 0.90*screenRatio/portraitRatio);
    width = portraitRatio*height/screenRatio;
    window = figure('Name', caption, 'NumberTitle', 'off', 'Color', 'w', ...
        'WindowStyle', 'normal', 'Units', 'normalized', ...
        'Position', [(1-width)/2 (1-height)/2 width height]);
    layout = tiledlayout(window, rows, columns, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(layout, caption);
end

function formatAxes(ax, ordinate)
    grid(ax, 'on');
    xlabel(ax, 't, с');
    ylabel(ax, ordinate);
    ax.FontSize = 9;
end

function label = deltaLabel(coefficients)
    terms = {};
    if coefficients(1) ~= 0
        terms{end+1} = sprintf('%g·δ(t)', coefficients(1));
    end
    if coefficients(2) ~= 0
        terms{end+1} = sprintf('%g·δ′(t)', coefficients(2));
    end
    label = strjoin(terms, ' + ');
end
