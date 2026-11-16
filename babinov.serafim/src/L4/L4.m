clc;
clear;
close all;

s = tf('s');

W_open = 1 / (0.05*s^4 + 0.1*s^3 + s^2 + s + 1);
W_closed = feedback(W_open, 1);

% Критерий Гурвица
H_open = [
    0.1  1    0  0
    0.05 1    1  0
    0    0.1  1  0
    0    0.05 1  1
];

H_closed = [
    0.1  1    0  0
    0.05 1    2  0
    0    0.1  1  0
    0    0.05 1  2
];

minors_open = [H_open(1, 1), det(H_open(1:2, 1:2)), ...
    det(H_open(1:3, 1:3)), det(H_open)];
minors_closed = [H_closed(1, 1), det(H_closed(1:2, 1:2)), ...
    det(H_closed(1:3, 1:3)), det(H_closed)];

disp('Матрица Гурвица разомкнутой системы:');
disp(H_open);
disp('Миноры разомкнутой системы:');
disp(minors_open);

disp('Матрица Гурвица замкнутой системы:');
disp(H_closed);
disp('Миноры замкнутой системы:');
disp(minors_closed);

disp('Полюса разомкнутой системы:');
disp(pole(W_open));
disp('Полюса замкнутой системы:');
disp(pole(W_closed));

% Критерий Найквиста
critical_x_limits = [-3 2];
critical_y_limits = [-2 2];

figure('Color', 'w', 'Theme', 'light');

subplot(1, 2, 1);
nyquist(W_open);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Полный годограф');

subplot(1, 2, 2);
nyquist(W_open);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim(critical_x_limits);
ylim(critical_y_limits);
title('Область критической точки');

sgtitle('Годограф Найквиста разомкнутой системы');

% Критерий Михайлова
w = linspace(0, 20, 5000);
D = polyval([0.05 0.1 1 1 2], 1i*w);
zoom_max_frequency = 5;

figure('Color', 'w', 'Theme', 'light');

subplot(1, 2, 1);
plot(real(D), imag(D), 'b', 'LineWidth', 1.5);
hold on;
plot(real(D(1)), imag(D(1)), 'go', 'MarkerFaceColor', 'g');
grid on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Полная кривая');

subplot(1, 2, 2);
initial_part = w <= zoom_max_frequency;
plot(real(D(initial_part)), imag(D(initial_part)), 'b', 'LineWidth', 1.5);
hold on;
plot(real(D(1)), imag(D(1)), 'go', 'MarkerFaceColor', 'g');
grid on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Начальный участок');

sgtitle('Кривая Михайлова замкнутой системы');

% Переходные процессы
t = 0:0.01:20;

figure('Color', 'w', 'Theme', 'light');
step(W_open, W_closed, t);
grid on;
legend('Разомкнутая', 'Замкнутая', 'Location', 'best');
xlabel('t, с');
ylabel('y(t)');
title('Переходные процессы');
