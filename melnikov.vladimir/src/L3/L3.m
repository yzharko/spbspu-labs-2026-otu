clear;
clc;
close all;

K0 = 1.2;
T = 0.3;
T1 = 0.5;
T2 = 0.6;
T1a = 1.5;

K_values = [0.6 1.2 2.4 3.6];
T_values = [0.15 0.3 0.6 1.2];

w = logspace(-2, 2, 500);
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
    'Форсирующее', ...
    'Изодромное'};

n = numel(names);

W = {1/s, s, tf(1), 1/(T*s+1), ...
    1/(T2^2*s^2+T1a*s+1), ...
    1/(T2^2*s^2+T1*s+1), ...
    1/(T*s^2+1), ...
    1/(s*(T*s+1)), ...
    s/(T*s+1), ...
    T*s+1, ...
    (T*s+1)/s};

W0 = cell(n, 1);

for i = 1:n
    W0{i} = K0 * W{i};
end

plotGroup(W0, 1:n, {}, w, 'Исходные характеристики', names);

WK = cell(n, length(K_values));
labelsK = cell(1, length(K_values));

for j = 1:length(K_values)
    K = K_values(j);
    labelsK{j} = sprintf('K = %.2g', K);

    for i = 1:n
        WK{i,j} = K * W{i};
    end
end

plotGroup(WK, 1:n, labelsK, w, 'Варьирование K', names);

idsT = [4 7 8 9 10 11];

WT = cell(length(idsT), length(T_values));
labelsT = cell(1, length(T_values));

for j = 1:length(T_values)
    Tv = T_values(j);
    labelsT{j} = sprintf('T = %.2g', Tv);

    for i = 1:length(idsT)

        switch idsT(i)
            case 4
                WT{i,j} = K0/(Tv*s+1);

            case 7
                WT{i,j} = K0/(Tv*s^2+1);

            case 8
                WT{i,j} = K0/(s*(Tv*s+1));

            case 9
                WT{i,j} = K0*s/(Tv*s+1);

            case 10
                WT{i,j} = K0*(Tv*s+1);

            case 11
                WT{i,j} = K0*(Tv*s+1)/s;
        end
    end
end

plotGroup(WT, idsT, labelsT, w, 'Варьирование T', names);

T1a_values = [1.2 1.5 1.8 2.4];
T2a_values = [0.3 0.4 0.5 0.6];

W5 = cell(2, length(T1a_values));

labels5_T1 = cell(1, length(T1a_values));
labels5_T2 = cell(1, length(T2a_values));

for j = 1:length(T1a_values)

    T1v = T1a_values(j);
    T2v = T2a_values(j);

    W5{1,j} = K0/(T2^2*s^2+T1v*s+1);

    W5{2,j} = K0/(T2v^2*s^2+T1a*s+1);

    labels5_T1{j} = sprintf('T1 = %.2g', T1v);
    labels5_T2{j} = sprintf('T2 = %.2g', T2v);
end

plotSecondOrder(W5, labels5_T1, labels5_T2, w, 5);

T1b_values = [0.2 0.35 0.5 0.8];
T2b_values = [0.35 0.5 0.6 0.9];

W6 = cell(2, length(T1b_values));

labels6_T1 = cell(1, length(T1b_values));
labels6_T2 = cell(1, length(T2b_values));

for j = 1:length(T1b_values)

    T1v = T1b_values(j);
    T2v = T2b_values(j);

    W6{1,j} = K0/(T2^2*s^2+T1v*s+1);

    W6{2,j} = K0/(T2v^2*s^2+T1*s+1);

    labels6_T1{j} = sprintf('T1 = %.2g', T1v);
    labels6_T2{j} = sprintf('T2 = %.2g', T2v);
end

plotSecondOrder(W6, labels6_T1, labels6_T2, w, 6);

function plotGroup(systems, ids, labels, w, heading, names)

    kinds = {'ЛАЧХ', 'ЛФЧХ', 'АФЧХ'};

    nSystems = size(systems, 1);
    nCurves = size(systems, 2);

    cols = 4;
    if nSystems <= 6
        cols = 3;
    end

    rows = ceil(nSystems / cols);
    colors = lines(nCurves);

    for kind = 1:3

        figure('Name', [heading ' - ' kinds{kind}], ...
            'NumberTitle', 'off', ...
            'Position', [80 80 1400 800]);

        tl = tiledlayout(rows, cols, ...
            'TileSpacing', 'compact', ...
            'Padding', 'compact');

        title(tl, [heading ' - ' kinds{kind}]);

        for i = 1:nSystems

            ax = nexttile(tl);
            hold(ax, 'on');

            h = gobjects(1, nCurves);

            for j = 1:nCurves
                h(j) = drawResponse(ax, systems{i,j}, ...
                    w, kind, colors(j,:), ids(i));
            end
            
            if kind <= 2
                set(ax, 'XScale', 'log');
                xlim(ax, [w(1) w(end)]);
            end

            grid(ax, 'on');
            title(ax, names{ids(i)}, ...
            'FontSize', 9);

            if nCurves > 1
                legend(ax, h, labels, ...
                    'Location', 'best', 'FontSize', 7);
            end

            hold(ax, 'off');
        end
    end
end

function plotSecondOrder(systems, labels1, labels2, w, number)

    kinds = {'ЛАЧХ', 'ЛФЧХ', 'АФЧХ'};
    paramNames = {'Изменение T1', 'Изменение T2'};
    labelSets = {labels1, labels2};

    nCurves = size(systems, 2);
    colors = lines(nCurves);

    figure('Name', sprintf('Звено %d - T1 и T2', number), ...
        'NumberTitle', 'off', ...
        'Position', [80 80 1400 800]);

    tl = tiledlayout(2, 3, ...
        'TileSpacing', 'compact', ...
        'Padding', 'compact');

    title(tl, sprintf('Звено %d: влияние T1 и T2', number));

    for i = 1:2
        for kind = 1:3

            ax = nexttile(tl, (i-1)*3 + kind);
            hold(ax, 'on');

            h = gobjects(1, nCurves);

            for j = 1:nCurves
                h(j) = drawResponse(ax, systems{i,j}, ...
                    w, kind, colors(j,:), number);
            end
            if kind <= 2
                set(ax, 'XScale', 'log');
                xlim(ax, [w(1) w(end)]);
            end
            grid(ax, 'on');

            title(ax, [paramNames{i} ' - ' kinds{kind}]);

            legend(ax, h, labelSets{i}, ...
                'Location', 'best', 'FontSize', 8);

            hold(ax, 'off');
        end
    end
end

function h = drawResponse(ax, sys, w, kind, color, id)

    % Консервативное звено имеет разрыв на резонансе
    if id == 7
        omega0 = max(abs(imag(pole(sys))));
        masks = {w < omega0, w > omega0};
    else
        masks = {true(size(w))};
    end

    if kind == 1 || kind == 2

        [mag, phase] = bode(sys, w);

        mag = squeeze(mag);
        phase = squeeze(phase);

        if kind == 1
            y = 20*log10(mag);
            ylabel(ax, 'L, дБ');
        else
            y = phase;
            ylabel(ax, 'Фаза, град');
        end

        y = y(:)';

        for k = 1:length(masks)
            mask = masks{k};

            if k == 1
                h = semilogx(ax, w(mask), y(mask), ...
                    'Color', color, 'LineWidth', 1.3);
            else
                semilogx(ax, w(mask), y(mask), ...
                    'Color', color, 'LineWidth', 1.3, ...
                    'HandleVisibility', 'off');
            end
        end

        xlabel(ax, 'Частота, рад/с');

    else

        [re, im] = nyquist(sys, w);

        re = squeeze(re);
        im = squeeze(im);

        re = re(:)';
        im = im(:)';
        
        if ismember(id, [1 2])
            re(:) = 0;
        elseif id == 7
            im(:) = 0;
        elseif ismember(id, [10 11])
            re(:) = mean(re);
        end

        % Усилительное звено отображается точкой
        if max(abs(re-re(1))) < 1e-10 && ...
                max(abs(im-im(1))) < 1e-10

            h = plot(ax, re(1), im(1), 'o', ...
                'Color', color, ...
                'MarkerFaceColor', color, ...
                'MarkerSize', 6);

        else
            for k = 1:length(masks)
                mask = masks{k};

                % Положительные частоты
                if k == 1
                    h = plot(ax, re(mask), im(mask), ...
                        'Color', color, 'LineWidth', 1.3);
                else
                    plot(ax, re(mask), im(mask), ...
                        'Color', color, 'LineWidth', 1.3, ...
                        'HandleVisibility', 'off');
                end

                plot(ax, re(mask), -im(mask), '--', ...
                    'Color', color, 'LineWidth', 1, ...
                    'HandleVisibility', 'off');
            end
        end

        xlabel(ax, 'Re W(j\omega)');
        ylabel(ax, 'Im W(j\omega)');
    end
end
