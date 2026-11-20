clc;
clear;
close all;

s = tf('s');
W_open = 10 / (2*s^3 + 2*s^2 + s + 1);
W_closed = 10 / (2*s^3 + 2*s^2 + s + 11);

% Построение годографа Найквиста
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

% Построение кривой Михайлова
[num, den] = tfdata(W_closed, 'v');
D = round(den, 10);

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

% Прямая проверка устойчивости

t = 0:0.01:10;
figure;
hold on;
grid on;
step(W_open, t);
step(W_closed, t);
title('Переходные процессы');
legend('Разомкнутая система', 'Замкнутая система');
xlabel('Время t, с');
ylabel('Выходной сигнал');
hold off;

% Полюса
syms s;
disp('Полюса разомкнутой системы'); 
eq = 2*s^3 + 2*s^2 + s + 1 == 0;
sol = solve(eq, s);
disp(double(sol));
disp('Полюса замкнутой системы');
eq = 2*s^3 + 2*s^2 + s + 11 == 0;
sol = solve(eq, s);
disp(double(sol));

