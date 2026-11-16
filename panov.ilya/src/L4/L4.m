
clc;
clear;
close all;

s = tf('s');
W_open = 10/(2*s^5 + 3*s^4 + 3*s^3 + 0.5*s^2 + 0.5*s + 1);
W_closed = 10/(2*s^5 + 3*s^4 + 3*s^3 + 0.5*s^2 + 0.5*s + 11);

H_open = [
    3 0.5 1 0 0
    2 3 0.5 0 0
    0 3 0.5 1 0
    0 2 3 0.5 0
    0 0 3 0.5 1
];

H_closed = [
    3 0.5 11 0 0
    2 3 0.5 0 0
    0 3 0.5 11 0
    0 2 3 0.5 0
    0 0 3 0.5 11
];

Delta_open = zeros(1, 5);
Delta_closed = zeros(1, 5);

fprintf('\nОпределители Гурвица:\n');
for k = 1:5
    Delta_open(k) = det(H_open(1:k, 1:k));
    Delta_closed(k) = det(H_closed(1:k, 1:k));
    fprintf('Delta_%d: разомкнутая = %.2f; замкнутая = %.2f\n', k, Delta_open(k), Delta_closed(k));
end

figure('Color', 'w', 'Theme', 'light');
nyquist(W_open);
grid on;
hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
title('Годограф Найквиста разомкнутой системы W(s)', 'Color', 'k');

figure('Color', 'w', 'Theme', 'light');
nyquist(W_open);
grid on;
hold on;
plot(-1, 0, 'r+', 'MarkerSize', 14, 'LineWidth', 2);
xlim([-3 2]);
ylim([-2 2]);
title('Область критической точки', 'Color', 'k');

[num, den] = tfdata(W_closed, 'v');
D = round(den, 10);
w = 0:0.02:20;
Djw = polyval(D, 1i*w);

figure('Color', 'w', 'Theme', 'light');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8);
hold on;
plot(real(Djw(1)), imag(Djw(1)), 'o', 'Color', [0 0.4 0], 'MarkerFaceColor', [0 0.4 0], 'MarkerSize', 8);
plot(real(Djw(end)), imag(Djw(end)), 's', 'Color', [0.7 0 0], 'MarkerFaceColor', [0.7 0 0], 'MarkerSize', 8);
text(real(Djw(1)), imag(Djw(1)), '  \omega = 0', 'Color', 'k');
text(real(Djw(end)), imag(Djw(end)), '  \omega = max', 'Color', 'k');
plot([0 0], [min(imag(Djw)) max(imag(Djw))], 'k:', 'LineWidth', 0.8);
plot([min(real(Djw)) max(real(Djw))], [0 0], 'k:', 'LineWidth', 0.8);
grid on;
box on;
xlabel('Re D(j\omega)', 'Color', 'k');
ylabel('Im D(j\omega)', 'Color', 'k');
title('Кривая Михайлова D(j\omega)', 'Color', 'k');
axis equal;

figure('Color', 'w', 'Theme', 'light');
plot(real(Djw), imag(Djw), 'b-', 'LineWidth', 1.8);
hold on;
plot(real(Djw(1)), imag(Djw(1)), 'o', 'Color', [0 0.4 0], 'MarkerFaceColor', [0 0.4 0], 'MarkerSize', 8);
plot([10 17], [0 0], 'k:', 'LineWidth', 0.8);
grid on;
box on;
xlabel('Re D(j\omega)', 'Color', 'k');
ylabel('Im D(j\omega)', 'Color', 'k');
title('Начальный участок кривой Михайлова', 'Color', 'k');
axis equal;
xlim([10 17]);
ylim([-1 2]);

p_open = pole(W_open);
p_closed = pole(W_closed);

fprintf('\nПолюсы разомкнутой системы:\n');
disp(p_open);
fprintf('Число полюсов с положительной вещественной частью: %d\n', sum(real(p_open) > 0));
fprintf('\nПолюсы замкнутой системы:\n');
disp(p_closed);
fprintf('Число полюсов с положительной вещественной частью: %d\n', sum(real(p_closed) > 0));
fprintf('\nПроверка isstable: 1 - устойчива, 0 - неустойчива.\n');
fprintf('Разомкнутая система: %d\n', isstable(W_open));
fprintf('Замкнутая система: %d\n', isstable(W_closed));

t = 0:0.01:20;

figure('Color', 'w', 'Theme', 'light');
subplot(2, 1, 1);
step(W_open, t);
grid on;
title('Переходный процесс разомкнутой системы', 'Color', 'k');

subplot(2, 1, 2);
step(W_closed, t);
grid on;
title('Переходный процесс замкнутой системы', 'Color', 'k');
