clc;
clear;
close all;

s = tf('s');
W_open = 1 / (0.05*s^4 + 0.1*s^3 + s^2 + s + 1);
W_closed = 1 / (0.05*s^4 + 0.1*s^3 + s^2 + s + 2);

steadyOpen = 1;
steadyClosed = 0.5;

wCross = sqrt(10);
[re, im] = nyquist(W_open);
re = [steadyOpen; re(:); 0];
im = [0; im(:); 0];

limits = [-1.2 1.2 -1.25 1.25; -3 2 -2 2];
titles = {'Годограф Найквиста разомкнутой системы W(s)', ...
    'Область критической точки'};
for i = 1:2
    figure('Color', 'w');
    plot(re, im, 'b-', re, -im, 'b--');
    hold on;
    plot(-1, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
    grid on;
    axis(limits(i, :));
    xlabel('Re W(jω)');
    ylabel('Im W(jω)');
    title(titles{i});
    legend('ω \geq 0', 'ω \leq 0', 'Критическая точка', ...
        'Location', 'best');
end

[num, den] = tfdata(W_closed, 'v');
D = den / num(end);
wAxis = [sqrt(10 - sqrt(60)), wCross, sqrt(10 + sqrt(60))];
w = unique([0:0.02:5, wAxis]);
Djw = polyval(D, 1i*w);

figure('Color', 'w');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.5);
hold on;
plot(real(Djw(1)), imag(Djw(1)), 'go', 'MarkerFaceColor', 'g');
plot(real(Djw(end)), imag(Djw(end)), 'rs', 'MarkerFaceColor', 'r');
text(real(Djw(1)), imag(Djw(1)), '  ω = 0');
text(real(Djw(end)), imag(Djw(end)), 'ω = 5  ', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');
xline(0, 'k:');
yline(0, 'k:');
grid on;
axis equal;
xlabel('Re D(jω)');
ylabel('Im D(jω)');
title('Кривая Михайлова D(jω)');

t = (0:0.001:15)';
yOpen = step(W_open, t);
yClosed = step(W_closed, t);

figure('Color', 'w');
plot(t, yOpen, 'b-', t, yClosed, 'r-', 'LineWidth', 1.4);
hold on;
yline(steadyOpen, 'b--', 'HandleVisibility', 'off');
yline(steadyClosed, 'r--', 'HandleVisibility', 'off');
grid on;
axis([0 15 0 1.25]);
xlabel('t, с');
ylabel('y(t)');
title('Переходные характеристики при единичном скачке');
legend('Разомкнутая САУ', 'Замкнутая САУ', 'Location', 'best');

disp('Полюсы разомкнутой системы');
disp(pole(W_open));
disp('Полюсы замкнутой системы');
disp(pole(W_closed));
[peakOpen, idxOpen] = max(yOpen);
[peakClosed, idxClosed] = max(yClosed);
fprintf('Разомкнутая САУ: предел %.3f, максимум %.4f при t = %.3f с.\n', ...
    steadyOpen, peakOpen, t(idxOpen));
fprintf('Замкнутая САУ: предел %.3f, максимум %.4f при t = %.3f с.\n', ...
    steadyClosed, peakClosed, t(idxClosed));
