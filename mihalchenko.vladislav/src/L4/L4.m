clc;
clear;
close all;

s = tf('s');
W_open = 1 / (8 * s^3 + 4 * s^2 + 2 * s + 1);
W_closed = feedback(W_open, 1);

[num_open, D_open] = tfdata(W_open, 'v');
[num_closed, D_closed] = tfdata(W_closed, 'v');
D_open = D_open / num_open(end);
D_closed = D_closed / num_closed(end);

out = fullfile(fileparts(mfilename('fullpath')), 'results_L4');
if ~isfolder(out)
    mkdir(out);
end

%% Критерий Гурвица и полюса
H_open = [D_open(2) D_open(4) 0;
    D_open(1) D_open(3) 0;
    0 D_open(2) D_open(4)];
H_closed = [D_closed(2) D_closed(4) 0;
    D_closed(1) D_closed(3) 0;
    0 D_closed(2) D_closed(4)];

minors_open = [H_open(1, 1), det(H_open(1:2, 1:2)), det(H_open)];
minors_closed = [H_closed(1, 1), det(H_closed(1:2, 1:2)), det(H_closed)];

disp('Полином разомкнутой системы:');
disp(D_open);
disp('Полином замкнутой системы:');
disp(D_closed);
disp('Матрица и миноры Гурвица разомкнутой системы:');
disp(H_open);
disp(minors_open);
disp('Матрица и миноры Гурвица замкнутой системы:');
disp(H_closed);
disp(minors_closed);
disp('Полюса разомкнутой системы:');
disp(pole(W_open));
disp('Полюса замкнутой системы:');
disp(pole(W_closed));

%% Годограф Найквиста разомкнутой системы
poleFrequency = 0.5;
frequencyGap = 0.001;
wLow = linspace(0, poleFrequency - frequencyGap, 1500);
wHigh = logspace(log10(poleFrequency + frequencyGap), 2, 2500);
nyquistLimits = [-6 6 -6 6; -1.3 0.6 -1.5 1.5];
nyquistTitles = {'Годограф Найквиста', 'Область критической точки'};
nyquistFiles = {'nyquist_overview.png', 'nyquist_zoom.png'};

for view = 1:2
    fig = figure('Visible', 'off', 'Theme', 'light', 'Color', 'w', ...
        'Position', [100 100 850 540]);
    nyquist(W_open, 'b', wLow);
    hold on;
    nyquist(W_open, 'b', wHigh);
    axis equal;
    axis(nyquistLimits(view, :));
    grid off;
    xlabel('Re W(j\omega)');
    ylabel('Im W(j\omega)');
    title(nyquistTitles{view});
    exportgraphics(fig, fullfile(out, nyquistFiles{view}), 'Resolution', 185);
    close(fig);
end

%% Кривая Михайлова замкнутой системы
w = unique([0:0.002:2, 0.5, 1 / sqrt(2)]);
Djw = polyval(D_closed, 1i * w);
axisFrequencies = [0, 0.5, 1 / sqrt(2)];
axisPoints = polyval(D_closed, 1i * axisFrequencies);
mikhailovLimits = [-15 3 -62 3; -0.6 2.2 -2 0.6];
mikhailovTitles = {'Кривая Михайлова', 'Начальный участок'};

fig = figure('Visible', 'off', 'Theme', 'light', 'Color', 'w', ...
    'Position', [100 100 1100 400]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for view = 1:2
    nexttile;
    plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.4);
    hold on;
    plot(real(axisPoints), imag(axisPoints), 'ro', 'MarkerFaceColor', 'r');
    xline(0, ':');
    yline(0, ':');
    axis(mikhailovLimits(view, :));
    grid on;
    xlabel('Re D(j\omega)');
    ylabel('Im D(j\omega)');
    title(mikhailovTitles{view});
end
exportgraphics(fig, fullfile(out, 'mikhailov.png'), 'Resolution', 185);
close(fig);

%% Переходные процессы при единичном скачке
time = (0:0.02:50)';

fig = figure('Visible', 'off', 'Theme', 'light', 'Color', 'w', ...
    'Position', [100 100 1100 400]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
step(W_open, 'b', time);
grid on;
xlabel('t');
ylabel('h(t)');
title('Разомкнутая система');
nexttile;
step(W_closed, 'r', time);
grid on;
xlabel('t');
ylabel('h(t)');
title('Замкнутая система');
exportgraphics(fig, fullfile(out, 'step.png'), 'Resolution', 185);
close(fig);

fprintf('Графики сохранены в %s\n', out);
