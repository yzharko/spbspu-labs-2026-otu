clear;
close all;
clc;

s = tf('s');
W_open = 1 / (0.1*s^3 + 0.1*s^2 + s + 1);
W_closed = feedback(W_open, 1);

figure('Color', 'w');
nyquist(W_open);
grid on;
hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-4 3]);
ylim([-4 4]);
title('Годограф Найквиста разомкнутой системы');
xlabel('Re W(j\omega)');
ylabel('Im W(j\omega)');

figure('Color', 'w');
nyquist(W_open);
grid on;
hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Область критической точки');
xlim([-1.5 0.5]);
ylim([-1 1]);
xlabel('Re W(j\omega)');
ylabel('Im W(j\omega)');

[~, den_closed] = tfdata(W_closed, 'v');
D = round(den_closed, 10);
w = 0:0.01:10;
Djw = polyval(D, 1i*w);

figure('Color', 'w');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8);
hold on;
plot(real(Djw(1)), imag(Djw(1)), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(real(Djw(end)), imag(Djw(end)), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
text(real(Djw(1)), imag(Djw(1)), '  \omega = 0', 'Color', 'g');
text(real(Djw(end)), imag(Djw(end)), '  \omega = 10', 'Color', 'r');
plot([0 0], [min(imag(Djw)) max(imag(Djw))], 'k:', 'LineWidth', 0.8);
plot([min(real(Djw)) max(real(Djw))], [0 0], 'k:', 'LineWidth', 0.8);
grid on;
box on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Кривая Михайлова замкнутой системы');
xlim([-10 3]);
ylim([-100 10]);

t = 0:0.01:20;

figure('Color', 'w');
subplot(1, 2, 1);
step(W_open, t);
grid on;
title('Разомкнутая система');
xlabel('t, с');
ylabel('y(t)');

subplot(1, 2, 2);
step(W_closed, t);
grid on;
title('Замкнутая система');
xlabel('t, с');
ylabel('y(t)');
