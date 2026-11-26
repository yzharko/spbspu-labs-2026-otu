clc; clear; close all;

s = tf('s');
W_open = 1 / (s^3 + 0.1*s^2 + 0.1*s + 1);
W_closed = feedback(W_open, 1);

[~, den_open] = tfdata(W_open, 'v');
[~, den_closed] = tfdata(W_closed, 'v');
A = round(den_open, 10);
D = round(den_closed, 10);

H_open = [A(2) A(4) 0; A(1) A(3) 0; 0 A(2) A(4)];
H_closed = [D(2) D(4) 0; D(1) D(3) 0; 0 D(2) D(4)];

fprintf('Миноры Гурвица (разомкнутая): D1 = %.2f, D2 = %.2f, D3 = %.2f\n', ...
    det(H_open(1,1)), det(H_open(1:2,1:2)), det(H_open));
fprintf('Миноры Гурвица (замкнутая):   D1 = %.2f, D2 = %.2f, D3 = %.2f\n\n', ...
    det(H_closed(1,1)), det(H_closed(1:2,1:2)), det(H_closed));

figure('Color', 'w');
nyquist(W_open);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Годограф Найквиста разомкнутой системы W(s)');

figure('Color', 'w');
nyquist(W_open);
grid on; hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-2.5 2]); ylim([-1.5 1.5]);
title('Область критической точки (-1; j0)');

w_max = 5.5;
w_step = 0.005;
w = 0:w_step:w_max;
Djw = polyval(D, 1i*w);

figure('Color', 'w');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8); hold on;
plot(real(Djw(1)), imag(Djw(1)), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot(real(Djw(end)), imag(Djw(end)), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

text(real(Djw(1)), imag(Djw(1)), '  \omega = 0', 'Color', [0 0.6 0]);
text(real(Djw(end)), imag(Djw(end)), sprintf('  \\omega = %.1f', w_max), 'Color', 'r');

xline(0, 'k:', 'LineWidth', 0.8);
yline(0, 'k:', 'LineWidth', 0.8);

grid on; box on;
xlabel('Re D(j\omega)');
ylabel('Im D(j\omega)');
title('Кривая Михайлова замкнутой системы D(j\omega)');

disp('Полюса разомкнутой системы:');
disp(pole(W_open));
disp('Полюса замкнутой системы:');
disp(pole(W_closed));

t_sim = 0:0.01:30;
figure('Color', 'w');
subplot(2, 1, 1);
step(W_open, t_sim);
grid on;
title('Переходный процесс разомкнутой системы');
ylabel('Амплитуда');

subplot(2, 1, 2);
step(W_closed, t_sim);
grid on;
title('Переходный процесс замкнутой системы');
ylabel('Амплитуда');
xlabel('Время, с');
