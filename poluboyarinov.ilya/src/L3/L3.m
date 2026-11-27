clear;
clc;
close all;

K = 1.8;
T = 0.6;
T1 = 0.9;
T2 = 1.2;
T11 = 3.0;

K_values = [0.6 1.2 1.8 3.6];
T_values = [0.3 0.6 0.9 1.2];
T1_values = [0.45 0.9 1.35 1.8];
T2_values = [0.6 0.9 1.2 1.5];
T2_aper_values = [0.6 0.9 1.2 1.4];
T11_values = [2.4 3.0 3.6 4.8];

frequency_min = 1e-2;
frequency_max = 1e2;
frequency_points = 1000;
w = logspace(log10(frequency_min), log10(frequency_max), frequency_points);
s = tf('s');

names = {'Интегрирующее', ...
    'Идеальное дифференцирующее', ...
    'Усилительное', ...
    'Апериодическое 1-го порядка', ...
    'Апериодическое 2-го порядка', ...
    'Колебательное', ...
    'Консервативное', ...
    'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', ...
    'Форсирующее 1-го порядка', ...
    'Изодромное'};

W = cell(11, 1);
W{1} = K/s;
W{2} = K*s;
W{3} = tf(K);
W{4} = K/(T*s+1);
W{5} = K/(T2^2*s^2+T11*s+1);
W{6} = K/(T2^2*s^2+T1*s+1);
W{7} = K/(T*s^2+1);
W{8} = K/(s*(T*s+1));
W{9} = K*s/(T*s+1);
W{10} = K*(T*s+1);
W{11} = K*(T*s+1)/s;

plotGroup(W, 1:11, {}, w, 'Исходные характеристики', names);

WK = cell(11, numel(K_values));
labelsK = cell(1, numel(K_values));

for j = 1:numel(K_values)
    labelsK{j} = sprintf('K = %.2g', K_values(j));

    for i = 1:11
        WK{i,j} = (K_values(j)/K)*W{i};
    end
end

plotGroup(WK, 1:11, labelsK, w, 'Влияние коэффициента K', names);

idsT = [4 7 8 9 10 11];
WT = cell(numel(idsT), numel(T_values));
labelsT = cell(1, numel(T_values));

for j = 1:numel(T_values)
    Tv = T_values(j);
    labelsT{j} = sprintf('T = %.2g', Tv);

    WT{1,j} = K/(Tv*s+1);
    WT{2,j} = K/(Tv*s^2+1);
    WT{3,j} = K/(s*(Tv*s+1));
    WT{4,j} = K*s/(Tv*s+1);
    WT{5,j} = K*(Tv*s+1);
    WT{6,j} = K*(Tv*s+1)/s;
end

plotGroup(WT, idsT, labelsT, w, 'Влияние постоянной времени T', names);

W2 = cell(4, numel(T2_aper_values));
labels2 = cell(1, numel(T2_aper_values));

for j = 1:numel(T2_aper_values)
    T2v = T2_aper_values(j);
    W2{1,j} = K/(T2v^2*s^2+T11*s+1);
    labels2{j} = sprintf('T2 = %.2g', T2v);
end

W2{2,1} = K/(T2^2*s^2+T11_values(1)*s+1);
W2{2,2} = K/(T2^2*s^2+T11_values(2)*s+1);
W2{2,3} = K/(T2^2*s^2+T11_values(3)*s+1);
W2{2,4} = K/(T2^2*s^2+T11_values(4)*s+1);
labels11 = arrayfun(@(x) sprintf('T11 = %.2g', x), ...
    T11_values, 'UniformOutput', false);

W2{3,1} = K/(T2^2*s^2+T1_values(1)*s+1);
W2{3,2} = K/(T2^2*s^2+T1_values(2)*s+1);
W2{3,3} = K/(T2^2*s^2+T1_values(3)*s+1);
W2{3,4} = K/(T2^2*s^2+T1_values(4)*s+1);
labels1 = arrayfun(@(x) sprintf('T1 = %.2g', x), ...
    T1_values, 'UniformOutput', false);

W2{4,1} = K/(T2_values(1)^2*s^2+T1*s+1);
W2{4,2} = K/(T2_values(2)^2*s^2+T1*s+1);
W2{4,3} = K/(T2_values(3)^2*s^2+T1*s+1);
W2{4,4} = K/(T2_values(4)^2*s^2+T1*s+1);
labels2_osc = arrayfun(@(x) sprintf('T2 = %.2g', x), ...
    T2_values, 'UniformOutput', false);

labelsSecondOrder = {labels2, labels11, labels1, labels2_osc};
idsSecondOrder = [5 5 6 6];
titlesSecondOrder = {'Апериодическое 2-го порядка: изменение T2', ...
    'Апериодическое 2-го порядка: изменение T11', ...
    'Колебательное звено: изменение T1', ...
    'Колебательное звено: изменение T2'};

plotSecondOrder(W2, idsSecondOrder, labelsSecondOrder, ...
    titlesSecondOrder, w);

function plotGroup(systems, ids, labels, w, heading, names)

    kinds = {'ЛАЧХ', 'ЛФЧХ', 'АФЧХ'};
    nSystems = size(systems, 1);
    nCurves = size(systems, 2);
    colors = lines(nCurves);
    nCols = 4;
    nRows = ceil(nSystems/nCols);

    for kind = 1:3
        figure('Name', [heading ' - ' kinds{kind}], ...
            'NumberTitle', 'off', ...
            'Position', [80 80 1400 850]);

        tl = tiledlayout(nRows, nCols, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        title(tl, [heading ' - ' kinds{kind}]);

        for i = 1:nSystems
            ax = nexttile(tl);
            hold(ax, 'on');
            h = gobjects(1, nCurves);

            for j = 1:nCurves
                h(j) = drawResponse(ax, systems{i,j}, w, ...
                    kind, colors(j,:), ids(i));
            end

            grid(ax, 'on');
            title(ax, names{ids(i)}, 'FontSize', 9);

            if kind <= 2
                set(ax, 'XScale', 'log');
                xlim(ax, [w(1) w(end)]);
                xlabel(ax, 'Частота, рад/с');
            else
                xlabel(ax, 'Re W(j\omega)');
                ylabel(ax, 'Im W(j\omega)');
            end

            if kind == 1
                ylabel(ax, 'L, дБ');
            elseif kind == 2
                ylabel(ax, 'Фаза, град');
            end

            if nCurves > 1 && ~isempty(labels)
                legend(ax, h, labels, ...
                    'Location', 'best', 'FontSize', 7);
            end

            hold(ax, 'off');
        end
    end
end

function plotSecondOrder(systems, ids, labels, headings, w)

    kinds = {'ЛАЧХ', 'ЛФЧХ', 'АФЧХ'};
    colors = lines(size(systems, 2));

    for kind = 1:3
        figure('Name', ['Звенья второго порядка - ' kinds{kind}], ...
            'NumberTitle', 'off', ...
            'Position', [80 80 1200 850]);

        tl = tiledlayout(2, 2, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        title(tl, ['Влияние постоянных времени - ' kinds{kind}]);

        for i = 1:4
            ax = nexttile(tl);
            hold(ax, 'on');
            h = gobjects(1, size(systems, 2));

            for j = 1:size(systems, 2)
                h(j) = drawResponse(ax, systems{i,j}, w, ...
                    kind, colors(j,:), ids(i));
            end

            grid(ax, 'on');
            title(ax, headings{i}, 'FontSize', 9);

            if kind <= 2
                set(ax, 'XScale', 'log');
                xlim(ax, [w(1) w(end)]);
                xlabel(ax, 'Частота, рад/с');
            else
                xlabel(ax, 'Re W(j\omega)');
                ylabel(ax, 'Im W(j\omega)');
            end

            if kind == 1
                ylabel(ax, 'L, дБ');
            elseif kind == 2
                ylabel(ax, 'Фаза, град');
            end

            legend(ax, h, labels{i}, ...
                'Location', 'best', 'FontSize', 8);
            hold(ax, 'off');
        end
    end
end

function h = drawResponse(ax, sys, w, kind, color, id)

    if id == 7
        p = pole(sys);
        omega0 = max(abs(imag(p)));
        relativeGap = 0.01;
        masks = {w < omega0*(1-relativeGap), ...
            w > omega0*(1+relativeGap)};
    else
        masks = {true(size(w))};
    end

    h = gobjects(1);

    for k = 1:numel(masks)
        currentW = w(masks{k});

        if isempty(currentW)
            continue;
        end

        if kind <= 2
            [mag, phase] = bode(sys, currentW);
            mag = squeeze(mag);
            phase = squeeze(phase);

            if kind == 1
                y = 20*log10(mag);
            else
                phase = phase(:)'*pi/180;
                y = unwrap(phase)*180/pi;
            end

            currentHandle = semilogx(ax, currentW, y(:)', ...
                'Color', color, 'LineWidth', 1.3);
        else
            [re, im] = nyquist(sys, currentW);
            re = squeeze(re);
            im = squeeze(im);
            re = re(:)';
            im = im(:)';

            currentHandle = plot(ax, re, im, ...
                'Color', color, 'LineWidth', 1.3);

            plot(ax, re, -im, '--', ...
                'Color', color, 'LineWidth', 1, ...
                'HandleVisibility', 'off');
        end

        if ~isgraphics(h)
            h = currentHandle;
        end
    end
end
