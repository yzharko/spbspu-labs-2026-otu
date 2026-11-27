clc;
clear;
close all;

s = tf('s');
W = 100 / (5*s^4 + 0.1*s^3 + 2*s^2 + 2*s + 1);
W_closed = feedback(W, 1);

[~, den_open]   = tfdata(W, 'v');
[~, den_closed] = tfdata(W_closed, 'v');
disp(roots(den_open));
disp(roots(den_closed));

figure('Theme', 'light', 'Color', 'w');
nyquist(W);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Годограф Найквиста разомкнутой системы W(s)');

figure('Theme', 'light', 'Color', 'w');
nyquist(W); grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-3 3]); ylim([-3 3]);
title('Область критической точки');

D = round(den_closed, 10);

w = 0:0.005:6;
Djw = polyval(D, 1i*w);

figure('Theme', 'light', 'Color', 'w');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8); hold on;
plot(real(Djw(1)), imag(Djw(1)), 'go', 'MarkerFaceColor', 'g', ...
     'MarkerSize', 8);
plot(real(Djw(end)), imag(Djw(end)), 'rs', 'MarkerFaceColor', 'r', ...
     'MarkerSize', 8);

text(real(Djw(1)),   imag(Djw(1)),   '  \omega = 0',   'Color', 'g');
text(real(Djw(end)), imag(Djw(end)), '  \omega = max', 'Color', 'r');

plot([0 0], [min(imag(Djw)) max(imag(Djw))], 'k:', 'LineWidth', 0.8);
plot([min(real(Djw)) max(real(Djw))], [0 0], 'k:', 'LineWidth', 0.8);

grid on; box on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Кривая Михайлова D(j\omega)');

t = 0:0.005:10;

figure('Theme', 'light', 'Color', 'w');
subplot(1, 2, 1);
step(W, t); grid on;
title('Разомкнутая система');

subplot(1, 2, 2);
step(W_closed, t); grid on;
title('Замкнутая система');
