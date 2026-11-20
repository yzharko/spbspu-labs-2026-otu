clc;
clear;
close all;

try, st = settings; st.matlab.appearance.figure.GraphicsTheme.TemporaryValue = "light"; catch, end
set(groot, 'defaultFigureColor', 'w', 'defaultAxesColor', 'w', ...
    'defaultAxesXColor', 'k', 'defaultAxesYColor', 'k', ...
    'defaultAxesGridColor', 'k', 'defaultTextColor', 'k');

s = tf('s');
W_open = 1 / (s^3 + 0.1*s^2 + 0.1*s + 1);
W_closed = feedback(W_open, 1);

figure('Color', 'w');
nyquist(W_open);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Годограф Найквиста разомкнутой системы W(s)');

figure('Color', 'w');
nyquist(W_open); grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-3 2]); ylim([-2 2]);
title('Область критической точки');

[num, den] = tfdata(W_closed, 'v');
D = round(den, 10);

w = 0:0.001:3;
Djw = polyval(D, 1i*w);

figure('Color', 'w');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8); hold on;
plot(real(Djw(1)), imag(Djw(1)), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(real(Djw(end)), imag(Djw(end)), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
text(real(Djw(1)), imag(Djw(1)), '  \omega = 0',   'Color', 'g');
text(real(Djw(end)), imag(Djw(end)), '  \omega = max', 'Color', 'r');
plot([0 0], [min(imag(Djw)) max(imag(Djw))], 'k:', 'LineWidth', 0.8);
plot([min(real(Djw)) max(real(Djw))], [0 0], 'k:', 'LineWidth', 0.8);
grid on; box on;
xlabel('Re D(j\omega)'); ylabel('Im D(j\omega)');
title('Кривая Михайлова D(j\omega)');

disp('Полюса замкнутой системы:'); disp(pole(W_closed));
disp('Полюса разомкнутой системы:'); disp(pole(W_open));

t = 0:0.01:40;
figure('Color', 'w');
step(W_open, W_closed, t); grid on;
legend('разомкнутая', 'замкнутая', 'Location', 'northwest');
title('Сравнение переходных процессов');
