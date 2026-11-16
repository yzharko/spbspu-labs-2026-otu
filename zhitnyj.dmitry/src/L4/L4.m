clc;
clear;
close all;

W_open = tf(1, [0.05, 0.1, 1, 1, 1]);
W_closed = feedback(W_open, 1);

[~, denOpen] = tfdata(W_open, 'v');
[~, denClosed] = tfdata(W_closed, 'v');
[hurwitzOpen, deltaOpen] = hurwitzCriteria(denOpen);
[hurwitzClosed, deltaClosed] = hurwitzCriteria(denClosed);
rootsOpen = pole(W_open);
rootsClosed = pole(W_closed);

disp('Коэффициенты характеристического полинома замкнутой системы');
disp(denClosed);
disp('Матрица Гурвица разомкнутой системы');
disp(hurwitzOpen);
disp('Матрица Гурвица замкнутой системы');
disp(hurwitzClosed);
disp(table((1:4)', deltaOpen, deltaClosed, ...
    'VariableNames', {'Определитель', 'Разомкнутая', 'Замкнутая'}));
disp(table(rootsOpen, rootsClosed, ...
    'VariableNames', {'Корни разомкнутой', 'Корни замкнутой'}));

omega = sort([0, logspace(-3, 3, 4001), sqrt(10)]);
nyquistFigure = figure('Name', 'Критерий Найквиста', ...
    'NumberTitle', 'off', 'Theme', 'light', 'Color', 'w');
layout = tiledlayout(nyquistFigure, 1, 2);

nexttile(layout, 1);
nyquist(W_open, omega);
nyquistChart = gca;
nyquistChart.Title.String = 'Годограф разомкнутой системы';
nyquistChart.XLabel.String = 'Re W(j\omega)';
nyquistChart.YLabel.String = 'Im W(j\omega)';
grid off;

nexttile(layout, 2);
nyquist(W_open, omega);
nyquistChart = gca;
nyquistChart.Title.String = 'Область критической точки';
nyquistChart.XLabel.String = 'Re W(j\omega)';
nyquistChart.YLabel.String = 'Im W(j\omega)';
nyquistChart.XLimits = [-1.5, 0.5];
nyquistChart.YLimits = [-1, 1];
grid off;

mikhailovFigure = figure('Name', 'Критерий Михайлова', ...
    'NumberTitle', 'off', 'Theme', 'light', 'Color', 'w');
layout = tiledlayout(mikhailovFigure, 1, 2);
frequencyLimits = [10, 5];
viewNames = ["Кривая Михайлова", "Начальный участок кривой"];
for view = 1:numel(frequencyLimits)
    w = linspace(0, frequencyLimits(view), 3001);
    characteristic = polyval(denClosed, 1i * w);
    ax = nexttile(layout);
    plot(ax, real(characteristic), imag(characteristic));
    hold(ax, 'on');
    plot(ax, real(characteristic(1)), imag(characteristic(1)), 'o');
    text(ax, real(characteristic(1)), imag(characteristic(1)), ...
        ' \omega = 0', 'Interpreter', 'tex', 'VerticalAlignment', 'bottom');
    xline(ax, 0, ':k');
    yline(ax, 0, ':k');
    grid(ax, 'on');
    axis(ax, 'equal');
    title(ax, viewNames(view));
    xlabel(ax, 'Re D(j\omega)', 'Interpreter', 'tex');
    ylabel(ax, 'Im D(j\omega)', 'Interpreter', 'tex');
end

time = 0:0.01:30;
[stepOpen, timeOpen] = step(W_open, time);
[stepClosed, timeClosed] = step(W_closed, time);
stepFigure = figure('Name', 'Переходные процессы', ...
    'NumberTitle', 'off', 'Theme', 'light', 'Color', 'w');
layout = tiledlayout(stepFigure, 1, 2);

ax = nexttile(layout);
plot(ax, timeOpen, stepOpen);
grid(ax, 'on');
title(ax, 'Разомкнутая система');
xlabel(ax, 't, с');
ylabel(ax, 'y(t)');

ax = nexttile(layout);
plot(ax, timeClosed, stepClosed);
grid(ax, 'on');
title(ax, 'Замкнутая система');
xlabel(ax, 't, с');
ylabel(ax, 'y(t)');

function [matrix, minors] = hurwitzCriteria(a)
    matrix = [
        a(2), a(4), 0,    0
        a(1), a(3), a(5), 0
        0,    a(2), a(4), 0
        0,    a(1), a(3), a(5)
    ];
    minors = zeros(4, 1);
    for order = 1:4
        minors(order) = det(matrix(1:order, 1:order));
    end
end
