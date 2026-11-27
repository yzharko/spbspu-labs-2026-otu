clc;
clear;
close all;

s = tf('s');

K0 = 4.0;
T0 = 1.5;
T1 = 1.6;
T2 = 2.0;
T1a = 5.0;

names = {'Интегрирующее', 'Дифференцирующее', 'Усилительное', ...
    'Апериодическое 1-го порядка', 'Апериодическое 2-го порядка', ...
    'Колебательное', 'Консервативное', 'Инерционное интегрирующее', ...
    'Реальное дифференцирующее', 'Форсирующее 1-го порядка', 'Изодромное'};

W = {K0/s, K0*s, tf(K0, 1), K0/(T0*s + 1), ...
    K0/(T2^2*s^2 + T1a*s + 1), K0/(T2^2*s^2 + T1*s + 1), ...
    K0/(T0*s^2 + 1), K0/(s*(T0*s + 1)), K0*s/(T0*s + 1), ...
    K0*(T0*s + 1), K0*(T0*s + 1)/s};

for i = 1:length(W)
    figure;
    bode(W{i});
    grid on;
    title(sprintf('№%d. %s. ЛАЧХ и ЛФЧХ', i, names{i}));

    figure;
    nyquist(W{i});
    grid on;
    title(sprintf('№%d. %s. АФЧХ', i, names{i}));
end

K_values = [1 4 8];
T_values = [0.1 0.5 1 2 5];
omega = logspace(-2, 2, 600);

WT = cell(length(W), length(T_values));
for i = 1:length(T_values)
    T = T_values(i);
    scale = T/T0;
    WT(:,i) = {W{1}; W{2}; W{3}; K0/(T*s + 1); ...
        K0/((scale*T2)^2*s^2 + scale*T1a*s + 1); ...
        K0/((scale*T2)^2*s^2 + scale*T1*s + 1); ...
        K0/(T*s^2 + 1); K0/(s*(T*s + 1)); K0*s/(T*s + 1); ...
        K0*(T*s + 1); K0*(T*s + 1)/s};
end

headings = {'ЛАЧХ при изменении K', 'ЛФЧХ при изменении K', ...
    'ЛАЧХ при изменении постоянных времени', 'ЛФЧХ при изменении постоянных времени'};

for panel = 1:4
    fig = figure('WindowStyle', 'normal', 'WindowState', 'maximized');
    theme(fig, 'light');
    layout = tiledlayout(3, 4, 'TileSpacing', 'compact');
    title(layout, headings{panel});
    xlabel(layout, 'Круговая частота, рад/с');
    isMagnitude = mod(panel, 2) == 1;

    if panel == 1
        values = K_values;
        labels = compose('K = %g', K_values);
    elseif panel == 2
        values = K0;
        labels = {'Фаза одинакова для всех K'};
    else
        values = T_values;
        labels = compose('T = %g с', T_values);
    end

    for n = 1:length(W)
        ax = nexttile;
        title(sprintf('Звено №%d', n));
        if panel >= 3 && n <= 3
            text(0.5, 0.5, 'Постоянных времени нет', 'HorizontalAlignment', 'center');
            axis off;
            continue;
        end

        hold(ax, 'on');
        ax.XScale = 'log';
        h = gobjects(1, length(values));
        for i = 1:length(values)
            if panel <= 2
                Wplot = W{n} * (values(i)/K0);
                T = T0;
            else
                Wplot = WT{n,i};
                T = values(i);
            end

            w = omega;
            valid = true(size(w));
            if n == 7
                poleFrequency = 1/sqrt(T);
                w = sort([w(w ~= poleFrequency), poleFrequency]);
                valid = w ~= poleFrequency;
            end
            [mag, phase] = bode(Wplot, w(valid));
            y = NaN(size(w));
            if isMagnitude
                y(valid) = 20*log10(mag(:));
            else
                y(valid) = phase(:);
            end

            h(i) = plot(ax, w, y, '-', 'Marker', 'none');
        end

        if isMagnitude
            ylabel('L, дБ');
        else
            ylabel('Фаза, °');
        end
        if panel == 1
            yline(0, ':');
            xline(1, ':');
        end
        xlim([omega(1) omega(end)]);
        grid on;
    end
    lgd = legend(h, labels, 'Orientation', 'horizontal');
    lgd.Layout.Tile = 'south';
end

