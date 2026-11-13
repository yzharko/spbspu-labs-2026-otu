clc;
clear;
close all;

s = tf('s');
W_opened = 1 / (8*s^3 + 4*s^2 + 2*s + 1);
W_closed = feedback(W_opened, 1);
display(W_closed);

%% Построение годографа Найквиста
figure('Color', 'w');
nyquist(W_opened);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Годограф Найквиста разомкнутой системы W(s)');
hold off;

%% Построение кривой Михайлова
[num, den] = tfdata(W_closed, 'v');
D = round(den, 10);
 
w = 0:0.01:1;              
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
t = 0:0.01:120;

[y_open, t_open] = step(W_opened, t);
[y_closed, t_closed] = step(W_closed, t);

figure('Color', 'w');
plot(t_open, y_open, 'b-', 'LineWidth', 1);
grid on;
title('Переходный процесс разомкнутой системы');
xlabel('Время, с'); ylabel('Амплитуда');

figure('Color', 'w');
plot(t_closed, y_closed, 'r-', 'LineWidth', 1);
grid on;
title('Переходный процесс замкнутой системы');
xlabel('Время, с'); ylabel('Амплитуда');

figure('Color', 'w');
plot(t_open, y_open, 'b-', 'LineWidth', 1); hold on;
plot(t_closed, y_closed, 'r--', 'LineWidth', 1);
grid on;
title('Сравнение переходных процессов');
xlabel('Время, с'); ylabel('Амплитуда');
legend('Разомкнутая система', 'Замкнутая система', 'Location', 'best');
xlim([0 50]);

[~, den_open] = tfdata(W_opened, 'v');
[~, den_closed] = tfdata(W_closed, 'v');

roots_open = roots(den_open);
roots_closed = roots(den_closed);

fprintf("Корни разомкнутой системы\n")
disp(roots_open);

fprintf("Корни замкнутой системы\n")
disp(roots_closed);
