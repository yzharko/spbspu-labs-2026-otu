clc;
clear;
close all;

Dopen = [1, 0.1, 0.1, 1];
Wopen = tf(1, Dopen);
Wclosed = feedback(Wopen, 1);
[~, Dclosed] = tfdata(Wclosed, 'v');

%% Критерий Гурвица и полюсы
disp('Определители Гурвица разомкнутой системы:');
disp(hurwitz_minors(Dopen));
disp('Определители Гурвица замкнутой системы:');
disp(hurwitz_minors(Dclosed));
disp('Полюсы разомкнутой системы:');
disp(pole(Wopen));
disp('Полюсы замкнутой системы:');
disp(pole(Wclosed));

%% Годограф Найквиста
w = [0, logspace(-3, 2, 1500)];
[re, im] = nyquist(Wopen, w);
figure('Name', 'Годограф Найквиста', 'NumberTitle', 'off', 'Theme', 'light', ...
    'Position', [70, 100, 900, 520]);
plot(re(:), im(:), 'b-', re(:), -im(:), 'r--', 'LineWidth', 1.5);
hold on;
plot(-1, 0, 'kx', 'MarkerSize', 10, 'LineWidth', 1.6);
xline(0, ':', 'HandleVisibility', 'off');
yline(0, ':', 'HandleVisibility', 'off');
grid on;
axis equal;
xlim([-1.15, 1.15]);
ylim([-0.68, 0.68]);
xlabel('Re W(j\omega)');
ylabel('Im W(j\omega)');
legend('\omega \geq 0', '\omega \leq 0', '(-1; 0)', 'Location', 'northwest');
title('Годограф Найквиста разомкнутой системы');

%% Кривая Михайлова
figure('Name', 'Кривая Михайлова', 'NumberTitle', 'off', 'Theme', 'light', ...
    'Position', [90, 110, 900, 440]);
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
wLimits = [6, 0.45];
for wMax = wLimits
    w = linspace(0, wMax, 1200);
    D = polyval(Dclosed, 1i*w);
    nexttile;
    plot(real(D), imag(D), 'b-', 'LineWidth', 1.5);
    hold on;
    plot(real(D(1)), imag(D(1)), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    xline(0, ':', 'HandleVisibility', 'off');
    yline(0, ':', 'HandleVisibility', 'off');
    grid on;
    xRange = [min(real(D)), max(real(D))];
    margin = 0.08 * diff(xRange);
    xlim(xRange + [-margin, margin]);
    xlabel('Re D(j\omega)');
    ylabel('Im D(j\omega)');
    title(sprintf('0 \\leq \\omega \\leq %g рад/с', wMax));
    legend('D(j\omega)', '\omega = 0', 'Location', 'northwest');
end

%% Ступенчатые отклики
t = 0:0.01:15;
figure('Name', 'Ступенчатые отклики', 'NumberTitle', 'off', 'Theme', 'light', ...
    'Position', [110, 120, 900, 530]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(t, step(Wopen, t), 'b-', 'LineWidth', 1.5);
grid on;
xlim([t(1), t(end)]);
ylabel('y(t)');
title('Разомкнутая система');
nexttile;
plot(t, step(Wclosed, t), 'r-', 'LineWidth', 1.5);
grid on;
xlim([t(1), t(end)]);
xlabel('Время, с');
ylabel('y(t)');
title('Замкнутая система');

function minors = hurwitz_minors(a)
    H = [a(2), a(4), 0;
         a(1), a(3), 0;
         0,    a(2), a(4)];
    minors = [H(1,1), det(H(1:2,1:2)), det(H)];
end
