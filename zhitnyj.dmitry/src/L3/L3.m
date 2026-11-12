clear; clc; close all;

p = struct('K', 1.5, 'T', 0.4, 'T1Ap', 2.0, ...
           'T1Osc', 0.6, 'T2', 0.8);
s = tf('s');
omega = logspace(-2, 2, 801).';

systems = createLinks(s, p);
names = [
    "Идеальное интегрирующее"
    "Идеальное дифференцирующее"
    "Усилительное"
    "Апериодическое первого порядка"
    "Апериодическое второго порядка"
    "Колебательное"
    "Консервативное"
    "Инерционное интегрирующее"
    "Реальное дифференцирующее"
    "Форсирующее первого порядка"
    "Изодромное"
];

for i = 1:numel(systems)
    fig = figure('Name', sprintf('%02d. %s', i, char(names(i))), ...
                 'NumberTitle', 'off', 'Theme', 'light', 'Color', 'w');
    layout = tiledlayout(fig, 1, 2);
    nexttile(layout, 1);
    plotArgs = {systems{i}, omega};
    if i == 7
        omega0 = 1 / sqrt(p.T);
        edges = omega0 * [0.995, 1.005];
        below = [omega(omega < edges(1)); edges(1)];
        above = [edges(2); omega(omega > edges(2))];
        plotArgs = {frd(systems{i}, below), frd(systems{i}, above)};
        phaseOptions = bodeoptions;
        phaseOptions.PhaseWrapping = 'on';
        phaseOptions.PhaseWrappingBranch = -180;
        bodeplot(plotArgs{:}, phaseOptions);
    else
        bode(systems{i}, omega);
    end
    bodeChart = gca;
    bodeChart.Title.String = 'ЛАЧХ и ЛФЧХ';
    bodeChart.YLabel.String = ["ЛАЧХ", "ЛФЧХ"];
    grid on;

    nexttile(layout, 2);
    if ismember(i, [1, 2, 10, 11])
        nyquistChart = nyquistplot(plotArgs{:});
        x = real(evalfr(systems{i}, 1i * omega(1)));
        nyquistChart.XLimits = x + abs(p.K) * [-1, 1];
    else
        nyquist(plotArgs{:});
        nyquistChart = gca;
    end
    nyquistChart.Title.String = 'АФЧХ';
    if i == 3
        nyquistChart.Responses(1).MarkerStyle = 'o';
        nyquistChart.XLimits = [-1, 2 * p.K + 1];
        nyquistChart.YLimits = [-1, 1];
    end
    if i == 7
        ylim(nyquistChart, [-1, 1]);
    end
    grid off;
end

Kvalues = p.K * [0.5, 1, 2];
for group = 1:2
    indices = (1 + 6 * (group - 1)):min(6 * group, numel(systems));
    figure('Name', 'Влияние коэффициента K', ...
           'NumberTitle', 'off', 'Theme', 'light', 'Color', 'w');
    layout = tiledlayout(2, 3, 'TileSpacing', 'compact');

    for i = indices
        ax = nexttile(layout);
        hold(ax, 'on');
        ax.XScale = 'log';
        for j = 1:numel(Kvalues)
            [mag, ~] = bode(systems{i} * (Kvalues(j) / p.K), omega);
            [w, L] = breakAtResonance(omega, 20 * log10(mag(:)), i, p.T);
            semilogx(ax, w, L, 'LineWidth', 1.4, ...
                'DisplayName', sprintf('K = %.2g', Kvalues(j)));
        end
        yline(ax, 0, ':k', 'HandleVisibility', 'off');
        grid(ax, 'on');
        title(ax, names(i));
        xlabel(ax, '\omega, рад/с', 'Interpreter', 'tex');
        ylabel(ax, 'L, дБ');
        if i == indices(1)
            legend(ax, 'show', 'Location', 'best');
        end
    end
end

cases = {
    4,  'T',     [0.2 0.4 0.8];
    7,  'T',     [0.2 0.4 0.8];
    8,  'T',     [0.2 0.4 0.8];
    9,  'T',     [0.2 0.4 0.8];
    10, 'T',     [0.2 0.4 0.8];
    11, 'T',     [0.2 0.4 0.8];
    5,  'T1Ap',  [1.8 2.0 2.4];
    5,  'T2',    [0.6 0.8 0.9];
    6,  'T1Osc', [0.4 0.6 0.9];
    6,  'T2',    [0.6 0.8 0.9];
};

for first = 1:2:size(cases, 1)
    figure('Name', 'Влияние постоянных времени', ...
           'NumberTitle', 'off', 'Theme', 'light', 'Color', 'w');
    layout = tiledlayout(2, 2, 'TileSpacing', 'compact');

    for row = 0:1
        k = first + row;
        i = cases{k, 1};
        field = cases{k, 2};
        values = cases{k, 3};

        axL = nexttile(layout);
        axP = nexttile(layout);
        hold(axL, 'on');
        hold(axP, 'on');
        set([axL, axP], 'XScale', 'log');

        for j = 1:numel(values)
            varied = p;
            varied.(field) = values(j);
            changedSystems = createLinks(s, varied);
            [mag, phase] = bode(changedSystems{i}, omega);
            [w, curves] = breakAtResonance(omega, ...
                [20 * log10(mag(:)), phase(:)], i, varied.T);

            lineName = sprintf('%s = %.2g', field, values(j));
            semilogx(axL, w, curves(:, 1), 'LineWidth', 1.4, ...
                     'DisplayName', lineName);
            semilogx(axP, w, curves(:, 2), 'LineWidth', 1.4);
        end

        grid(axL, 'on');
        grid(axP, 'on');
        title(axL, names(i) + ': ЛАЧХ');
        title(axP, names(i) + ': ЛФЧХ');
        xlabel(axL, '\omega, рад/с', 'Interpreter', 'tex');
        xlabel(axP, '\omega, рад/с', 'Interpreter', 'tex');
        ylabel(axL, 'L, дБ');
        ylabel(axP, 'Фаза, ^{\circ}', 'Interpreter', 'tex');
        legend(axL, 'show', 'Location', 'best');
    end
end

function systems = createLinks(s, p)
    systems = {
        p.K / s
        p.K * s
        tf(p.K)
        p.K / (p.T * s + 1)
        p.K / (p.T2^2 * s^2 + p.T1Ap * s + 1)
        p.K / (p.T2^2 * s^2 + p.T1Osc * s + 1)
        p.K / (p.T * s^2 + 1)
        p.K / (s * (p.T * s + 1))
        p.K * s / (p.T * s + 1)
        p.K * (p.T * s + 1)
        p.K * (p.T * s + 1) / s
    };
end

function [w, y] = breakAtResonance(w, y, linkIndex, T)
    if linkIndex == 7
        split = find(w < 1 / sqrt(T), 1, 'last');
        if ~isempty(split) && split < numel(w)
            w = [w(1:split); NaN; w(split+1:end)];
            y = [y(1:split, :); NaN(1, size(y, 2)); y(split+1:end, :)];
        end
    end
end
