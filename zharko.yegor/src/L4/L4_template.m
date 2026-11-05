clc;
clear;
close all;

s = tf('s');
W_closed = ((3*s + 5)*(s^2 + 3*s + 2)) / (s^3 + 4*s^2 + 75*s + 112);

%% Построение годографа Найквиста
figure('Color', 'w');
nyquist(W_closed);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Годограф Найквиста замкнутой системы W(s)');

figure('Color', 'w');
nyquist(W_closed); grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-3 2]); ylim([-2 2]);
title('Область критической точки');

%% Построение кривой Михайлова
[num, den] = tfdata(W_closed, 'v');
D = round(den, 10);
 
% Задаём диапазон частот (при необходимости увеличиваем верхнюю границу)
w = 0:0.02:20;              
Djw = polyval(D, 1i*w);
 
figure('Color', 'w');
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
axis equal;

%% Прямая проверка устойчивости
???
