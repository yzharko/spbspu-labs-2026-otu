clc;
clear;
close all;

s = tf('s');
W_open = 3 / (0.1*s^3 + 0.01*s^2 + 0.1*s + 1);
W_closed = feedback(W_open, 1);

[~, a] = tfdata(W_open, 'v');
[~, b] = tfdata(W_closed, 'v');

H_open = [a(2) a(4) 0; a(1) a(3) 0; 0 a(2) a(4)];
H_closed = [b(2) b(4) 0; b(1) b(3) 0; 0 b(2) b(4)];

disp('Полюса разомкнутой системы:');
disp(pole(W_open));
disp('Матрица Гурвица разомкнутой системы:');
disp(H_open);
disp('Определители Гурвица разомкнутой системы:');
disp([H_open(1,1), det(H_open(1:2, 1:2)), det(H_open)]);

disp('Полюса замкнутой системы:');
disp(pole(W_closed));
disp('Матрица Гурвица замкнутой системы:');
disp(H_closed);
disp('Определители Гурвица замкнутой системы:');
disp([H_closed(1,1), det(H_closed(1:2, 1:2)), det(H_closed)]);

figure('Color', 'w');
nyquist(W_open);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
title('Годограф Найквиста разомкнутой системы W(s)');

figure('Color', 'w');
nyquist(W_open);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
xlim([-3 2]); ylim([-1.5 1.5]);
title('Область критической точки');

w = 0:0.01:50;
D = polyval(b, 1i*w);

figure('Color', 'w');
plot(real(D), imag(D), 'b-', 'LineWidth', 1.8);
hold on; grid on; box on;
plot(real(D(1)), imag(D(1)), 'go', 'MarkerFaceColor', 'g');
plot(real(D(end)), imag(D(end)), 'rs', 'MarkerFaceColor', 'r');
xlabel('Re D(j\omega)'); ylabel('Im D(j\omega)');
title('Кривая Михайлова D(j\omega)');

w2 = 0:0.001:2;
D2 = polyval(b, 1i*w2);

figure('Color', 'w');
plot(real(D2), imag(D2), 'b-', 'LineWidth', 1.8);
hold on; grid on; box on;
plot(real(D2(1)), imag(D2(1)), 'go', 'MarkerFaceColor', 'g');
xlabel('Re D(j\omega)'); ylabel('Im D(j\omega)');
title('Начальный участок кривой Михайлова');

t = 0:0.005:5;
figure('Color', 'w');
subplot(2, 1, 1);
step(W_open, t);
grid on;
title('Переходный процесс разомкнутой системы');
xlabel('t, c'); ylabel('h(t)');

subplot(2, 1, 2);
step(W_closed, t);
grid on;
title('Переходный процесс замкнутой системы');
xlabel('t, c'); ylabel('h(t)');
