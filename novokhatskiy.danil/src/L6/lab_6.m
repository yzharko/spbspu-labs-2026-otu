clc;
clear;
close all;

%% 1. Объект и исходная замкнутая система
s = tf('s');
Go = 1 / ((s + 0.1) * (s + 0.5) * (s + 2));
tLimit = 2.0;
overshootLimit = 10;
settlingBand = 0.02;
G0 = feedback(Go, 1);

%% 2. Второй метод Зиглера — Николса
% D(s) = s^3 + a2*s^2 + a1*s + a0. Для P-регулятора
% строка s^1 таблицы Рауса обнуляется при a2*a1 = a0 + Kcrit.
[~, denominator] = tfdata(Go, 'v');
a2 = denominator(2);
a1 = denominator(3);
a0 = denominator(4);
Kcrit = a2 * a1 - a0;
wcrit = sqrt(a1);
Tcrit = 2 * pi / wcrit;
Gcrit = feedback(Kcrit * Go, 1);

Kp = 0.6 * Kcrit;
Ti = 0.5 * Tcrit;
Td = 0.125 * Tcrit;
Ki = Kp / Ti;
Kd = Kp * Td;
Czn = pid(Kp, Ki, Kd);

fprintf('Kкр = %.4f, ωкр = %.4f рад/с, Tкр = %.4f с\n', Kcrit, wcrit, Tcrit);
fprintf('Полюса системы при критическом усилении:\n');
disp(pole(Gcrit));
fprintf('Зиглер — Николс: Ti = %.4f с, Td = %.4f с\n', Ti, Td);

%% 3. Автонастройка и анализ трёх замкнутых систем
[Cauto, autoInfo] = pidtune(Go, 'PID');
Lzn = Czn * Go;
Lauto = Cauto * Go;
Gzn = feedback(Lzn, 1);
Gauto = feedback(Lauto, 1);
caseNames = ["Без регулятора"; "Зиглер — Николс"; "pidtune"];
loops = {Go; Lzn; Lauto};
systems = {G0; Gzn; Gauto};

fprintf('Коэффициенты ПИД в параллельной форме:\n');
gains = table(caseNames(2:3), [Czn.Kp; Cauto.Kp], ...
    [Czn.Ki; Cauto.Ki], [Czn.Kd; Cauto.Kd], ...
    'VariableNames', {'Method', 'Kp', 'Ki', 'Kd'});
disp(gains);
fprintf('Частота среза при автонастройке: %.4f рад/с\n\n', ...
    autoInfo.CrossoverFrequency);

settling = zeros(3, 1);
overshoot = zeros(3, 1);
finalValue = zeros(3, 1);
gainMarginDb = zeros(3, 1);
phaseMarginDeg = zeros(3, 1);
stable = false(3, 1);
for i = 1:numel(systems)
    sys = systems{i};
    metrics = stepinfo(sys, 'SettlingTimeThreshold', settlingBand);
    settling(i) = metrics.SettlingTime;
    overshoot(i) = metrics.Overshoot;
    finalValue(i) = dcgain(sys);
    stable(i) = isstable(sys);
    [gainMargin, phaseMarginDeg(i)] = margin(loops{i});
    gainMarginDb(i) = 20 * log10(gainMargin);
    fprintf('%s. Полюса замкнутой системы:\n', char(caseNames(i)));
    disp(pole(sys));
end

staticError = 1 - finalValue;
timeOK = settling <= tLimit;
overshootOK = overshoot <= overshootLimit;
results = table(caseNames, stable, finalValue, staticError, settling, ...
    overshoot, gainMarginDb, phaseMarginDeg, timeOK, overshootOK, ...
    'VariableNames', {'Method', 'Stable', 'FinalValue', 'StaticError', ...
    'Settling_s', 'Overshoot_pct', 'GainMargin_dB', 'PhaseMargin_deg', ...
    'TimeOK', 'OvershootOK'});
disp(results);

%% 4. Переходные характеристики при единичном скачке
blue = [0.10 0.34 0.58];
orange = [0.72 0.29 0.08];
grey = [0.68 0.70 0.73];
fontName = 'Arial';

fig0 = figure('Name', 'Исходная система', 'Theme', 'light', 'Color', 'w', ...
    'Position', [100 100 920 430]);
ax0 = axes(fig0, 'FontName', fontName, 'FontSize', 13);
hold(ax0, 'on');
t = linspace(0, 30, 2401);
y = step(G0, t);
band = finalValue(1) * [1 - settlingBand, 1 + settlingBand];
fill(ax0, [0 t(end) t(end) 0], band([1 1 2 2]), grey, ...
    'FaceAlpha', 0.24, 'EdgeColor', 'none', 'DisplayName', 'Полоса 2 %');
plot(ax0, t, y, 'Color', blue, 'LineWidth', 2.1, 'DisplayName', 'Без регулятора');
timeLabel = strrep(sprintf('t_п = %.2f с', settling(1)), '.', ',');
xline(ax0, settling(1), '--', timeLabel, ...
    'Interpreter', 'tex', 'LabelVerticalAlignment', 'middle', 'HandleVisibility', 'off');
yline(ax0, 1, ':', 'Задание', 'HandleVisibility', 'off');
xlim(ax0, [0 t(end)]);
ylim(ax0, [0 1.4]);
grid(ax0, 'on');
box(ax0, 'on');
xlabel(ax0, 'Время, с');
ylabel(ax0, 'Выход y(t)');
legend(ax0, 'Location', 'northeast');

figCrit = figure('Name', 'Критическое усиление', 'Theme', 'light', 'Color', 'w', ...
    'Position', [120 120 920 430]);
axCrit = axes(figCrit, 'FontName', fontName, 'FontSize', 13);
hold(axCrit, 'on');
t = linspace(0, 4*Tcrit, 2001);
y = step(Gcrit, t);
plot(axCrit, t, y, 'Color', blue, 'LineWidth', 2.1);
yline(axCrit, dcgain(Gcrit), ':');
grid(axCrit, 'on');
box(axCrit, 'on');
xlim(axCrit, [0 t(end)]);
ylim(axCrit, [0 2.03]);
xticks(axCrit, 0:2:22);
xlabel(axCrit, 'Время, с');
ylabel(axCrit, 'Выход y(t)');

figPID = figure('Name', 'Сравнение ПИД-регуляторов', 'Theme', 'light', 'Color', 'w', ...
    'Position', [140 140 920 430]);
axPID = axes(figPID, 'FontName', fontName, 'FontSize', 13);
hold(axPID, 'on');
tEnd = 5 * ceil(1.2 * max(settling(2:3)) / 5);
t = linspace(0, tEnd, 4001);
band = [1 - settlingBand, 1 + settlingBand];
hBand = fill(axPID, [0 tEnd tEnd 0], band([1 1 2 2]), grey, ...
    'FaceAlpha', 0.24, 'EdgeColor', 'none', 'DisplayName', 'Полоса 2 %');
hZN = plot(axPID, t, step(Gzn, t), 'Color', orange, 'LineWidth', 2.0, ...
    'DisplayName', 'Зиглер — Николс');
hAuto = plot(axPID, t, step(Gauto, t), '--', 'Color', blue, 'LineWidth', 2.2, ...
    'DisplayName', 'pidtune');
xline(axPID, tLimit, '--k', sprintf('t = %g с', tLimit), ...
    'LabelVerticalAlignment', 'middle', 'HandleVisibility', 'off');
yline(axPID, 1 + overshootLimit/100, ':k', ...
    sprintf('Предел %g %%', overshootLimit), ...
    'HandleVisibility', 'off');
grid(axPID, 'on');
box(axPID, 'on');
xlim(axPID, [0 tEnd]);
ylim(axPID, [0 1.65]);
xlabel(axPID, 'Время, с');
ylabel(axPID, 'Выход y(t)');
legend(axPID, [hZN hAuto hBand], 'Location', 'northeast');
