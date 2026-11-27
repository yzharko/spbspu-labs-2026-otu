clear;
clc;
close all;

numerator = 1;
denominatorOpen = [8 4 2 1];
denominatorClosed = denominatorOpen;
denominatorClosed(end) = denominatorClosed(end) + numerator;

openSystem = tf(numerator, denominatorOpen);
closedSystem = feedback(openSystem, 1, -1);
openPoles = pole(openSystem);
closedPoles = pole(closedSystem);
poleTolerance = 1e-8;

fprintf('D разомкнутой системы: %s\n', poly2str(denominatorOpen, 's'));
fprintf('D замкнутой системы: %s\n', poly2str(denominatorClosed, 's'));

[openHurwitz, openMinors] = hurwitzCubic(denominatorOpen);
[closedHurwitz, closedMinors] = hurwitzCubic(denominatorClosed);

fprintf('\nКритерий Гурвица\n');
fprintf('Матрица разомкнутой системы:\n');
disp(openHurwitz);
fprintf('Главные миноры: %g, %g, %g\n', openMinors);
fprintf('Полюса разомкнутой системы:\n');
disp(round(openPoles, 6));
fprintf('Разомкнутая система на границе устойчивости.\n\n');

fprintf('Матрица замкнутой системы:\n');
disp(closedHurwitz);
fprintf('Главные миноры: %g, %g, %g\n', closedMinors);
fprintf('Полюса замкнутой системы:\n');
disp(round(closedPoles, 6));
fprintf('Замкнутая система неустойчива.\n');

imaginaryPole = openPoles(abs(real(openPoles)) < poleTolerance & imag(openPoles) > poleTolerance);
poleFrequency = imag(imaginaryPole);
indentRadius = 0.01;
outerRadius = 30;

omegaBelow = linspace(0, poleFrequency - indentRadius, 1001).';
omegaAbove = logspace(log10(poleFrequency + indentRadius), log10(outerRadius), 2400).';
responseBelow = [dcgain(openSystem); nyquistValues(openSystem, omegaBelow(2:end))];
responseAbove = nyquistValues(openSystem, omegaAbove);

indentAngles = linspace(-pi/2, pi/2, 601).';
indentPoints = 1i*poleFrequency + indentRadius * exp(1i*indentAngles);
indentResponse = numerator ./ polyval(denominatorOpen, indentPoints);

outerAngles = linspace(pi/2, -pi/2, 401).';
outerPoints = outerRadius * exp(1i*outerAngles);
outerResponse = numerator ./ polyval(denominatorOpen, outerPoints);

positivePath = [responseBelow; indentResponse; responseAbove];
nyquistPath = [conj(flipud(positivePath)); positivePath(2:end); outerResponse];
phaseAroundCritical = unwrap(angle(1 + nyquistPath));
clockwiseEncirclements = round((phaseAroundCritical(1) - phaseAroundCritical(end)) / (2*pi));
openRightPoles = nnz(real(openPoles) > poleTolerance);
closedRightPoles = openRightPoles + clockwiseEncirclements;

assert(closedRightPoles == nnz(real(closedPoles) > poleTolerance), ...
    'Результат критерия Найквиста не совпал с подсчётом правых полюсов.');
fprintf('\nКритерий Найквиста\n');
fprintf('Правых полюсов разомкнутой системы: %d\n', openRightPoles);
fprintf('Охватов точки (-1, 0) по часовой стрелке: %d\n', clockwiseEncirclements);
fprintf('Правых полюсов замкнутой системы: %d\n', closedRightPoles);

positiveAxis = [responseBelow; NaN; responseAbove];
frequencyResponse = [conj(flipud(positiveAxis)); positiveAxis(2:end)];
arcResponse = [conj(flipud(indentResponse)); NaN; indentResponse; NaN; outerResponse];

figure('Name', 'Критерий Найквиста');
tiles = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tiles, 'Годограф Найквиста разомкнутой системы');

ax = nexttile(tiles);
plotNyquist(ax, frequencyResponse, arcResponse);
title(ax, 'Полный контур');
legend(ax, 'Location', 'best');

ax = nexttile(tiles);
plotNyquist(ax, frequencyResponse, arcResponse);
title(ax, 'Область точки (-1, 0)');
xlim(ax, [-2.2 1.8]);
ylim(ax, [-2 2]);

omegaMikhailov = linspace(0, 20, 4001).';
mikhailov = polyval(denominatorClosed, 1i*omegaMikhailov);
mikhailovPhase = unwrap(angle(mikhailov));
angleChange = rad2deg(mikhailovPhase(end) - mikhailovPhase(1));
requiredAngle = 90 * (numel(denominatorClosed) - 1);

fprintf('\nКритерий Михайлова\n');
fprintf('Изменение аргумента до %.0f рад/с: %.2f град. Для устойчивости: +%g град в пределе.\n', ...
    omegaMikhailov(end), angleChange, requiredAngle);
fprintf('Замкнутая система неустойчива.\n');

figure('Name', 'Критерий Михайлова');
ax = axes;
plot(ax, real(mikhailov), imag(mikhailov), 'LineWidth', 1.5);
hold(ax, 'on');

startPoint = mikhailov(1);
plot(ax, real(startPoint), imag(startPoint), 'o', ...
    'Color', [0.2 0.6 0.3], 'MarkerFaceColor', [0.2 0.6 0.3]);
text(ax, real(startPoint) - 0.08, imag(startPoint) - 0.08, '\omega = 0', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');

xline(ax, 0, ':');
yline(ax, 0, ':');
axis(ax, 'equal');
xlim(ax, [-2 2.25]);
ylim(ax, [-2 2]);
grid(ax, 'on');
xlabel(ax, 'Re D(j\omega)');
ylabel(ax, 'Im D(j\omega)');
title(ax, 'Кривая Михайлова замкнутой системы');

time = (0:0.05:70).';
[openStep, openTime] = step(openSystem, time);
[closedStep, closedTime] = step(closedSystem, time);

figure('Name', 'Переходные процессы');
tiles = tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tiles, 'Отклики на единичную ступеньку');

ax = nexttile(tiles);
plotStep(ax, openTime, openStep, 'Разомкнутая система');

ax = nexttile(tiles);
plotStep(ax, closedTime, closedStep, 'Замкнутая система');

function [matrix, minors] = hurwitzCubic(denominator)
    a0 = denominator(1);
    a1 = denominator(2);
    a2 = denominator(3);
    a3 = denominator(4);
    matrix = [a1 a3 0; a0 a2 0; 0 a1 a3];
    minors = [a1, det(matrix(1:2, 1:2)), det(matrix)];
end

function values = nyquistValues(system, frequencies)
    [realPart, imaginaryPart] = nyquist(system, frequencies);
    values = complex(realPart(:), imaginaryPart(:));
end

function plotNyquist(ax, frequencies, arcs)
    plot(ax, real(frequencies), imag(frequencies), 'LineWidth', 1.5, ...
        'DisplayName', 'W(j\omega)');
    hold(ax, 'on');
    plot(ax, real(arcs), imag(arcs), '--', 'LineWidth', 1.2, ...
        'DisplayName', 'Дуги контура');
    plot(ax, -1, 0, 'r+', 'MarkerSize', 10, 'LineWidth', 1.7, ...
        'DisplayName', 'Точка (-1, 0)');
    axis(ax, 'equal');
    grid(ax, 'on');
    xlabel(ax, 'Re W(s)');
    ylabel(ax, 'Im W(s)');
end

function plotStep(ax, time, values, caption)
    plot(ax, time, values, 'LineWidth', 1.5);
    grid(ax, 'on');
    xlabel(ax, 't, с');
    ylabel(ax, 'y(t)');
    title(ax, caption);
end
