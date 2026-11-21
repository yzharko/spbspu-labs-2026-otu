clc;
clear;
delete(findall(groot, 'Type', 'figure', 'Tag', 'L3FrequencyPlot'));
s = tf('s');

K0 = 1.0;
T0 = 0.2;
T1_0 = 0.4;
T2_0 = 0.5;
K_values = [0.5 1 2 3 5];
T_values = [0.1 T0 0.5 1 2 5];
omega = logspace(-2, 2, 1200);
showFigures = true;
outputDir = fullfile(fileparts(mfilename('fullpath')), 'L3_graphs');
if ~isfolder(outputDir)
    mkdir(outputDir);
end

T1_ap = max(T1_0, 2*T2_0);
secondOrder = @(K,T) K/(T(2)^2*s^2 + T(1)*s + 1);
links = {
    'Идеальное интегрирующее',       @(K,~) K/s,                 []
    'Идеальное дифференцирующее',    @(K,~) K*s,                 []
    'Усилительное',                 @(K,~) tf(K),              []
    'Апериодическое 1-го порядка',  @(K,T) K/(T*s + 1),        T0
    'Апериодическое 2-го порядка',  secondOrder,                [T1_ap T2_0]
    'Колебательное',                secondOrder,                [T1_0 T2_0]
    'Консервативное',               @(K,T) K/(T*s^2 + 1),      T0
    'Инерционное интегрирующее',    @(K,T) K/(s*(T*s + 1)),    T0
    'Реальное дифференцирующее',    @(K,T) K*s/(T*s + 1),      T0
    'Форсирующее 1-го порядка',     @(K,T) K*(T*s + 1),        T0
    'Изодромное',                   @(K,T) K*(T*s + 1)/s,      T0
};
T1_values = {[1 1.25 1.5 2 2.5], [0.1 0.25 T1_0 0.7 0.9]};
T2_values = {[0.2 0.3 0.4 0.45 T2_0], [0.25 0.35 T2_0 0.75 1]};

W = cell(1, size(links, 1));
poleFreq = zeros(size(W));
poleFreq(7) = 1/sqrt(T0);
for i = 1:numel(W)
    name = links{i, 1};
    model = links{i, 2};
    nominal = links{i, 3};
    W{i} = model(K0, nominal);
    plotResponse(name, W(i), {}, omega, poleFreq(i), i, fullfile(outputDir, sprintf('base_%02d.png', i)), showFigures);

    gains = K_values;
    if i == 6
        gains = sort([gains 0.9]);
    end
    WK = arrayfun(@(K) model(K, nominal), gains, 'UniformOutput', false);
    plotResponse([name ' — K'], WK, compose('K = %g', gains), omega, repmat(poleFreq(i), size(gains)), ...
        0, fullfile(outputDir, sprintf('K_%02d.png', i)), showFigures);

    if isempty(nominal)
        continue;
    end
    values = {T_values};
    parameters = {'T'};
    if numel(nominal) == 2
        values = {T1_values{i-4}, T2_values{i-4}};
        parameters = {'T1', 'T2'};
    end
    for j = 1:numel(values)
        varied = values{j};
        times = repmat(nominal, numel(varied), 1);
        times(:, j) = varied(:);
        WT = arrayfun(@(k) model(K0, times(k, :)), 1:numel(varied), 'UniformOutput', false);
        poles = (poleFreq(i) > 0)./sqrt(varied);
        parameter = parameters{j};
        plotResponse([name ' — ' parameter], WT, compose([parameter ' = %g'], varied), omega, poles, ...
            0, fullfile(outputDir, sprintf('%s_%02d.png', parameter, i)), showFigures);
    end
end

fprintf('\nT = %4.2f: 1/T = %6.3f рад/с, 1/sqrt(T) = %6.3f рад/с\n', [T_values; 1./T_values; 1./sqrt(T_values)]);
for i = 1:numel(W)
    fprintf('\n%2d. %s\n', i, links{i, 1});
    for K = K_values
        crossing = gainCrossings((K/K0)*W{i});
        result = strtrim(sprintf('%.3f  ', crossing));
        if i == 3 && K == 1
            result = 'все частоты';
        elseif isempty(crossing)
            result = 'нет';
        end
        fprintf('K = %3g: %s\n', K, result);
    end
end

function plotResponse(name, systems, labels, omega, poles, number, file, showFigures)
    sizePx = [725 545];
    box = [0.13 0.55 0.81 0.38];
    if number > 0
        sizePx = [1120 510];
        box = [0.08 0.57 0.40 0.34];
    end
    fig = figure('Name', name, 'NumberTitle', 'off', 'Tag', 'L3FrequencyPlot', ...
        'WindowStyle', 'normal', 'MenuBar', 'none', 'ToolBar', 'none', 'NextPlot', 'add', ...
        'Visible', 'on', 'Color', 'w', 'Units', 'pixels', 'Position', [60 60 sizePx]);
    magAx = plotAxes(fig, box, 'ЛАЧХ', '', 'L, дБ');
    phaseAx = plotAxes(fig, [box(1) 0.14 box(3) 0.31], 'ЛФЧХ', '\omega, рад/с', '\psi, град');
    set([magAx phaseAx], 'XScale', 'log', 'XMinorGrid', 'off', 'XLim', omega([1 end]));
    magAx.XTickLabel = [];
    yline(magAx, 0, ':', 'HandleVisibility', 'off');

    for j = 1:numel(systems)
        [w, valid] = frequencyGrid(omega, poles(j));
        [amplitude, degrees] = bode(systems{j}, w(valid));
        level = nan(size(w));
        phase = level;
        level(valid) = 20*log10(amplitude(:));
        phase(valid) = degrees(:);
        curve = plot(magAx, w, level, 'LineWidth', 1.4);
        if poles(j) > 0
            phase(valid) = -180*(w(valid) > poles(j));
            xline(phaseAx, poles(j), ':', 'Color', curve.Color);
        end
        plot(phaseAx, w, phase, 'Color', curve.Color, 'LineWidth', 1.4);
    end
    if number > 0
        ax = plotAxes(fig, [0.61 0.22 0.35 0.69], 'АФЧХ', 'Re W(j\omega)', 'Im W(j\omega)');
        drawNyquist(ax, systems{1}, omega, poles(1), number);
    else
        legend(magAx, labels, 'Location', 'northoutside', 'NumColumns', 3, 'TextColor', 'k', 'Box', 'off');
    end
    frame = getframe(fig);
    imwrite(frame.cdata, file);
    if ~showFigures
        delete(fig);
    end
end

function ax = plotAxes(fig, box, name, xText, yText)
    colors = [.12 .35 .61; .72 .25 .15; .15 .45 .32; .50 .30 .61; .61 .43 .07; .04 .52 .64];
    ax = axes(fig, 'Units', 'normalized', 'Position', box, 'Color', 'w', 'ColorOrder', colors, ...
        'XColor', 'k', 'YColor', 'k', 'Interactions', [], 'Toolbar', [], 'Box', 'on', 'NextPlot', 'add');
    disableDefaultInteractivity(ax);
    grid(ax, 'on');
    title(ax, name, 'Color', 'k');
    xlabel(ax, xText, 'Color', 'k');
    ylabel(ax, yText, 'Color', 'k');
end

function drawNyquist(ax, system, omega, pole, number)
    if ~ismember(number, [1 8 11])
        omega = [0 omega];
    end
    [w, valid] = frequencyGrid(omega, pole);
    [re, im] = nyquist(system, w(valid));
    u = nan(size(w));
    v = u;
    u(valid) = re(:);
    v(valid) = im(:);
    if number == 3
        plot(ax, u(1), v(1), 'b.', 'MarkerSize', 18);
    elseif number == 7
        v(valid) = 0;
        plot(ax, u, v, 'b', 'LineWidth', 1.4);
        axis(ax, [-6 6 -0.75 0.75]);
    else
        plot(ax, u, v, 'b', u, -v, 'r--', 'LineWidth', 1.4);
        legend(ax, {'\omega > 0', '\omega < 0'}, 'Location', 'southoutside', 'Orientation', 'horizontal', ...
            'TextColor', 'k', 'Box', 'off');
    end
    if ismember(number, [1 2 8 10 11])
        ylim(ax, [-5 5]);
    end
    axis(ax, 'square');
end

function [w, valid] = frequencyGrid(omega, pole)
    w = omega;
    if pole > 0
        w = sort([omega, pole*[0.99 0.999 1 1.001 1.01]]);
    end
    valid = pole == 0 | w ~= pole;
end

function omega0 = gainCrossings(system)
    [num, den] = tfdata(system, 'v');
    n = max(numel(num), numel(den));
    num = [zeros(1, n-numel(num)), num];
    den = [zeros(1, n-numel(den)), den];
    signs = (-1).^(n-1:-1:0);
    difference = conv(num, num.*signs) - conv(den, den.*signs);
    x = roots(difference(1:2:end).*signs);
    x = real(x(abs(imag(x)) < 1e-10 & real(x) > 0));
    omega0 = sort(sqrt(x));
end

